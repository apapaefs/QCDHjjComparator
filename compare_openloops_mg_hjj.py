#!/usr/bin/env python3
"""Compare OpenLoops and MadGraph ME2 values for gg -> Hjj subprocesses."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DEFAULT_SQRT_S = 13600.0
DEFAULT_MH = 125.0
DEFAULT_MJJ_MIN = 20.0
DEFAULT_PTJ = 20.0
DEFAULT_YMAX = 5.0
DEFAULT_MU = 125.0
DEFAULT_ALPHA_S = 0.11264802949303165

TSV_COLUMNS = [
    "event",
    "subprocess",
    "shat",
    "mg_raw",
    "ol_raw",
    "mg_me2",
    "ol_me2",
    "ol_over_mg",
    "mg_so_raw",
    "mg_retcode",
    "mg_acc",
    "ol_acc",
]


@dataclass(frozen=True)
class Momentum:
    pdg_id: int
    energy: float
    px: float
    py: float
    pz: float
    mass: float = 0.0

    def boost_z(self, y: float) -> "Momentum":
        c = math.cosh(y)
        s = math.sinh(y)
        return Momentum(
            self.pdg_id,
            c * self.energy + s * self.pz,
            self.px,
            self.py,
            s * self.energy + c * self.pz,
            self.mass,
        )

    def as_openloops(self) -> list[float]:
        return [self.energy, self.px, self.py, self.pz, self.mass]

    def as_madgraph(self) -> list[float]:
        return [self.energy, self.px, self.py, self.pz]


@dataclass(frozen=True)
class PhaseSpacePoint:
    event: int
    subprocess: str
    x1: float
    x2: float
    shat: float
    momenta: list[Momentum]


@dataclass(frozen=True)
class MGResult:
    raw: float
    so_raw: float
    retcode: int
    accuracy: float


@dataclass(frozen=True)
class OpenLoopsResult:
    raw: float
    accuracy: float


def lam(a: float, b: float, c: float) -> float:
    return max(0.0, a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


def direction(rng: random.Random) -> tuple[float, float, float]:
    c = 2.0 * rng.random() - 1.0
    s = math.sqrt(max(0.0, 1.0 - c * c))
    phi = 2.0 * math.pi * rng.random()
    return s * math.cos(phi), s * math.sin(phi), c


def boost(p: Momentum, beta: tuple[float, float, float]) -> Momentum:
    bx, by, bz = beta
    b2 = bx * bx + by * by + bz * bz
    if b2 == 0.0:
        return p
    gamma = 1.0 / math.sqrt(1.0 - b2)
    bp = bx * p.px + by * p.py + bz * p.pz
    fac = ((gamma - 1.0) * bp / b2) + gamma * p.energy
    return Momentum(
        p.pdg_id,
        gamma * (p.energy + bp),
        p.px + fac * bx,
        p.py + fac * by,
        p.pz + fac * bz,
        p.mass,
    )


def pt(p: Momentum) -> float:
    return math.hypot(p.px, p.py)


def rapidity(p: Momentum) -> float:
    num = p.energy + p.pz
    den = p.energy - p.pz
    if num <= 0.0 or den <= 0.0:
        return math.copysign(float("inf"), p.pz)
    return 0.5 * math.log(num / den)


def invariant_mass(momenta: list[Momentum]) -> float:
    e = sum(p.energy for p in momenta)
    px = sum(p.px for p in momenta)
    py = sum(p.py for p in momenta)
    pz = sum(p.pz for p in momenta)
    m2 = e * e - px * px - py * py - pz * pz
    # Boosted massless momenta can leave tiny roundoff-level residual masses.
    scale = max(e * e, px * px + py * py + pz * pz, 1.0)
    if abs(m2) < 1.0e-12 * scale:
        m2 = 0.0
    return math.sqrt(max(0.0, m2))


def passes_cuts(final_lab: list[Momentum], ptj: float, ymax: float, mjj_min: float) -> bool:
    jets = [p for p in final_lab if p.pdg_id != 25]
    if len(jets) != 2:
        return False
    return (
        min(pt(j) for j in jets) >= ptj
        and max(abs(rapidity(j)) for j in jets) < ymax
        and invariant_mass(jets) >= mjj_min
    )


def make_trial_point(
    event: int,
    subprocess: str,
    final_pdgs: list[int],
    rng: random.Random,
    sqrt_s: float,
    mh: float,
    mjj_min: float,
) -> tuple[PhaseSpacePoint, list[Momentum]]:
    if len(final_pdgs) != 2:
        raise ValueError("final_pdgs must contain exactly two final-state parton PDG ids")

    # Sample x1,x2 as tau/y, then build the partonic 2 -> 3 point by
    # sequentially decaying gg -> H + Q and Q -> parton + parton.
    s_had = sqrt_s * sqrt_s
    tau_min = (mh + mjj_min) ** 2 / s_had
    log_tau_min = math.log(tau_min)
    log_tau_width = -log_tau_min
    tau = math.exp(log_tau_min + rng.random() * log_tau_width)
    y_max = 0.5 * math.log(1.0 / tau)
    y = (2.0 * rng.random() - 1.0) * y_max
    x1 = math.sqrt(tau) * math.exp(y)
    x2 = math.sqrt(tau) * math.exp(-y)
    shat = tau * s_had
    ecm = math.sqrt(shat)

    s_pair_min = mjj_min * mjj_min
    s_pair_max = (ecm - mh) ** 2
    s_pair = s_pair_min + rng.random() * (s_pair_max - s_pair_min)
    m_pair = math.sqrt(s_pair)

    nx, ny, nz = direction(rng)
    pabs = math.sqrt(lam(shat, mh * mh, s_pair)) / (2.0 * ecm)
    higgs = Momentum(25, (shat + mh * mh - s_pair) / (2.0 * ecm), pabs * nx, pabs * ny, pabs * nz, mh)
    pair = Momentum(0, (shat - mh * mh + s_pair) / (2.0 * ecm), -pabs * nx, -pabs * ny, -pabs * nz, m_pair)

    px, py, pz = direction(rng)
    e_parton = 0.5 * m_pair
    p1_rf = Momentum(final_pdgs[0], e_parton, e_parton * px, e_parton * py, e_parton * pz, 0.0)
    p2_rf = Momentum(final_pdgs[1], e_parton, -e_parton * px, -e_parton * py, -e_parton * pz, 0.0)
    beta_pair = (pair.px / pair.energy, pair.py / pair.energy, pair.pz / pair.energy)
    p1 = boost(p1_rf, beta_pair)
    p2 = boost(p2_rf, beta_pair)

    # The provider-facing order is fixed to the current MG/OpenLoops convention:
    # incoming gluons first, then H, then the two selected final-state partons.
    initial = [
        Momentum(21, 0.5 * ecm, 0.0, 0.0, 0.5 * ecm, 0.0),
        Momentum(21, 0.5 * ecm, 0.0, 0.0, -0.5 * ecm, 0.0),
    ]
    final_cm = [higgs, p1, p2]
    point = PhaseSpacePoint(event, subprocess, x1, x2, shat, initial + final_cm)
    final_lab = [p.boost_z(y) for p in final_cm]
    return point, final_lab


def generate_accepted_points(
    subprocess: str,
    final_pdgs: list[int],
    n: int,
    seed: int,
    sqrt_s: float = DEFAULT_SQRT_S,
    mh: float = DEFAULT_MH,
    mjj_min: float = DEFAULT_MJJ_MIN,
    ptj: float = DEFAULT_PTJ,
    ymax: float = DEFAULT_YMAX,
    max_attempts: int | None = None,
) -> list[PhaseSpacePoint]:
    if n < 0:
        raise ValueError("n must be non-negative")
    if n == 0:
        return []
    attempts_limit = max_attempts if max_attempts is not None else max(1000, 1000 * n)
    rng = random.Random(seed)
    accepted: list[PhaseSpacePoint] = []
    attempts = 0
    while len(accepted) < n and attempts < attempts_limit:
        attempts += 1
        trial, final_lab = make_trial_point(
            len(accepted) + 1,
            subprocess,
            final_pdgs,
            rng,
            sqrt_s,
            mh,
            mjj_min,
        )
        if passes_cuts(final_lab, ptj, ymax, mjj_min):
            accepted.append(trial)
    if len(accepted) != n:
        raise RuntimeError(f"generated {len(accepted)} accepted points after {attempts} attempts; requested {n}")
    return accepted


def load_config(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def select_subprocess(config: dict[str, Any], name: str) -> dict[str, Any]:
    subprocesses = config.get("subprocesses", {})
    if name not in subprocesses:
        available = ", ".join(sorted(subprocesses)) or "<none>"
        raise ValueError(f"unknown subprocess '{name}'; available subprocesses: {available}")
    spec = dict(subprocesses[name])
    if "final_pdgs" not in spec:
        raise ValueError(f"subprocess '{name}' is missing final_pdgs")
    return spec


def write_mg_points(points: list[PhaseSpacePoint], path: Path, alpha_s: float, mu: float, category: str) -> None:
    with path.open("w") as out:
        for point in points:
            vals: list[Any] = [point.event, category, alpha_s, mu]
            for p in point.momenta:
                vals.extend(p.as_madgraph())
            out.write(" ".join(str(v) for v in vals) + "\n")


def read_mg_results(path: Path) -> dict[int, MGResult]:
    results: dict[int, MGResult] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            cols = line.split()
            if len(cols) < 6:
                raise ValueError(f"malformed MadGraph result line in {path}: {line.rstrip()}")
            event = int(cols[0])
            results[event] = MGResult(
                raw=float(cols[2]),
                so_raw=float(cols[3]),
                retcode=int(cols[4]),
                accuracy=float(cols[5]),
            )
    return results


def _read_existing(path: Path) -> bytes | None:
    return path.read_bytes() if path.exists() else None


def _restore_existing(path: Path, content: bytes | None) -> None:
    if content is None:
        if path.exists():
            path.unlink()
    else:
        path.write_bytes(content)


def validate_mg_paths(mg_dir: Path, executable: Path) -> None:
    if not mg_dir.exists():
        raise FileNotFoundError(
            f"Configured mg_dir does not exist: {mg_dir}. "
            "Do not run compare_openloops_mg_hjj.example.json directly; "
            "copy compare_openloops_mg_hjj.example.json to compare_openloops_mg_hjj.json "
            "and set mg_dir to the generated MadGraph SubProcesses/PV... directory."
        )
    if not mg_dir.is_dir():
        raise FileNotFoundError(f"Configured mg_dir is not a directory: {mg_dir}")
    if not executable.exists():
        raise FileNotFoundError(
            f"MadGraph evaluator executable not found: {executable}. "
            "Copy hgg_mg_eval.f into the generated subprocess directory and run "
            "'make -f Makefile -f /path/to/QCDHjjComparator/mg_eval.mk hgg_mg_eval'."
        )


def run_madgraph_evaluator(
    subprocess_config: dict[str, Any],
    points: list[PhaseSpacePoint],
    alpha_s: float,
    mu: float,
    category: str,
) -> dict[int, MGResult]:
    mg_dir = Path(subprocess_config["mg_dir"]).expanduser().resolve()
    executable = mg_dir / subprocess_config.get("mg_executable", "hgg_mg_eval")
    input_path = mg_dir / subprocess_config.get("mg_input", "hgg_mg_points.dat")
    output_path = mg_dir / subprocess_config.get("mg_output", "hgg_mg_eval.out")
    validate_mg_paths(mg_dir, executable)

    # Existing standalone MG drivers use fixed filenames; preserve any user's
    # current input/output files byte-for-byte around this temporary evaluation.
    old_input = _read_existing(input_path)
    old_output = _read_existing(output_path)
    try:
        write_mg_points(points, input_path, alpha_s, mu, category)
        if output_path.exists():
            output_path.unlink()
        proc = subprocess.run(
            [str(executable)],
            cwd=mg_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"MadGraph evaluator failed with {proc.returncode}\n{proc.stdout}")
        if not output_path.exists():
            raise RuntimeError(f"MadGraph evaluator did not write {output_path}\n{proc.stdout}")
        return read_mg_results(output_path)
    finally:
        _restore_existing(input_path, old_input)
        _restore_existing(output_path, old_output)


def evaluate_openloops(
    openloops_root: Path,
    subprocess_config: dict[str, Any],
    points: list[PhaseSpacePoint],
    alpha_s: float,
    mu: float,
    mh: float,
) -> dict[int, OpenLoopsResult]:
    root = openloops_root.expanduser().resolve()
    process_name = subprocess_config["openloops_process"]
    amptype = subprocess_config.get("openloops_amptype", "ls")

    old_cwd = Path.cwd()
    # The bundled OpenLoops Python wrapper expects to be imported from the
    # OpenLoops install root so its relative lib/proclib paths resolve.
    sys.path.insert(0, str(root / "pyol" / "tools"))
    try:
        os.chdir(root)
        import openloops  # type: ignore

        openloops.set_parameter("alpha_s", float(alpha_s))
        openloops.set_parameter("mu", float(mu))
        openloops.set_parameter("mass(25)", float(mh))
        proc = openloops.Process(process_name, amptype)
        results: dict[int, OpenLoopsResult] = {}
        for point in points:
            psp = [p.as_openloops() for p in point.momenta]
            me = proc.evaluate(psp)
            # For loop-induced Hjj channels the quantity of interest is the
            # loop-squared finite piece; fall back to loop.finite for other amps.
            raw = me.loop2.finite if hasattr(me, "loop2") else me.loop.finite
            results[point.event] = OpenLoopsResult(raw=float(raw), accuracy=float(me.acc))
        return results
    finally:
        os.chdir(old_cwd)


def make_comparison_rows(
    subprocess: str,
    subprocess_config: dict[str, Any],
    points: list[PhaseSpacePoint],
    mg_results: dict[int, MGResult],
    openloops_results: dict[int, OpenLoopsResult],
) -> list[dict[str, Any]]:
    mg_factor = float(subprocess_config.get("mg_factor", 1.0))
    ol_factor = float(subprocess_config.get("openloops_factor", 1.0))
    rows: list[dict[str, Any]] = []
    for point in points:
        if point.event not in mg_results:
            raise ValueError(f"missing MadGraph result for event {point.event}")
        if point.event not in openloops_results:
            raise ValueError(f"missing OpenLoops result for event {point.event}")
        mg = mg_results[point.event]
        ol = openloops_results[point.event]
        # Keep raw provider values visible, and put all convention corrections
        # behind explicit config factors for auditability.
        mg_me2 = mg.raw * mg_factor
        ol_me2 = ol.raw * ol_factor
        if not math.isfinite(mg_me2) or mg_me2 == 0.0:
            raise ValueError(f"invalid MadGraph ME2 for event {point.event}: {mg_me2}")
        if not math.isfinite(ol_me2):
            raise ValueError(f"invalid OpenLoops ME2 for event {point.event}: {ol_me2}")
        rows.append(
            {
                "event": point.event,
                "subprocess": subprocess,
                "shat": point.shat,
                "mg_raw": mg.raw,
                "ol_raw": ol.raw,
                "mg_me2": mg_me2,
                "ol_me2": ol_me2,
                "ol_over_mg": ol_me2 / mg_me2,
                "mg_so_raw": mg.so_raw,
                "mg_retcode": mg.retcode,
                "mg_acc": mg.accuracy,
                "ol_acc": ol.accuracy,
            }
        )
    return rows


def write_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: _format_tsv_value(row[col]) for col in TSV_COLUMNS})


def _format_tsv_value(value: Any) -> Any:
    if isinstance(value, float):
        return f"{value:.17e}"
    return value


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subprocess", required=True)
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--seed", type=int, default=260616)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--sqrt-s", type=float, default=DEFAULT_SQRT_S)
    ap.add_argument("--mh", type=float, default=DEFAULT_MH)
    ap.add_argument("--mjj-min", type=float, default=DEFAULT_MJJ_MIN)
    ap.add_argument("--ptj", type=float, default=DEFAULT_PTJ)
    ap.add_argument("--ymax", type=float, default=DEFAULT_YMAX)
    ap.add_argument("--mu", type=float, default=DEFAULT_MU)
    ap.add_argument("--alpha-s", type=float, default=DEFAULT_ALPHA_S)
    ap.add_argument("--max-attempts", type=int)
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config = load_config(args.config)
    subprocess_config = select_subprocess(config, args.subprocess)
    points = generate_accepted_points(
        subprocess=args.subprocess,
        final_pdgs=[int(p) for p in subprocess_config["final_pdgs"]],
        n=args.n,
        seed=args.seed,
        sqrt_s=args.sqrt_s,
        mh=args.mh,
        mjj_min=args.mjj_min,
        ptj=args.ptj,
        ymax=args.ymax,
        max_attempts=args.max_attempts,
    )
    mg_results = run_madgraph_evaluator(subprocess_config, points, args.alpha_s, args.mu, args.subprocess)
    openloops_results = evaluate_openloops(
        Path(config["openloops_root"]),
        subprocess_config,
        points,
        args.alpha_s,
        args.mu,
        args.mh,
    )
    rows = make_comparison_rows(args.subprocess, subprocess_config, points, mg_results, openloops_results)
    write_tsv(rows, args.output)
    print(f"wrote {args.output} with {len(rows)} comparison rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
