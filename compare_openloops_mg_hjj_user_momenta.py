#!/usr/bin/env python3
"""Compare OpenLoops and MadGraph for one user-supplied gg -> Hjj point.

The input point is deliberately plain text: five rows of E px py pz in the
provider order used by the existing MG helper, namely g, g, H, parton, parton.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path


DEFAULT_ALPHA_S = 0.11264802949303165
DEFAULT_MH = 125.0
DEFAULT_MU = 125.0

FIELDS = [
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


def read_momenta(path, mh=DEFAULT_MH):
    """Read one text momentum file and append the masses expected by OpenLoops."""
    rows = []
    for line in Path(path).read_text().splitlines():
        # Allow blank lines and comments so hand-written point files stay tidy.
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = [float(x) for x in line.split()]
        if len(parts) != 4:
            raise ValueError("each momentum line must contain exactly: E px py pz")
        rows.append(parts)
    if len(rows) != 5:
        raise ValueError(f"expected exactly 5 momentum lines, got {len(rows)}")
    masses = [0.0, 0.0, float(mh), 0.0, 0.0]
    return [row + [mass] for row, mass in zip(rows, masses)]


def shat(momenta):
    """Compute the partonic invariant mass squared from the two incoming legs."""
    e = momenta[0][0] + momenta[1][0]
    px = momenta[0][1] + momenta[1][1]
    py = momenta[0][2] + momenta[1][2]
    pz = momenta[0][3] + momenta[1][3]
    return e * e - px * px - py * py - pz * pz


def load_config(path):
    """Load the JSON file that maps subprocess names to provider settings."""
    with Path(path).open() as handle:
        return json.load(handle)


def get_subprocess(config, name):
    """Return the config block for one named subprocess, with a useful error."""
    try:
        return dict(config["subprocesses"][name])
    except KeyError as exc:
        names = ", ".join(sorted(config.get("subprocesses", {}))) or "<none>"
        raise ValueError(f"unknown subprocess '{name}'; available: {names}") from exc


def write_mg_input(path, momenta, alpha_s, mu, subprocess_name):
    """Write the one-line input format used by the simple MG evaluator."""
    values = [1, subprocess_name, alpha_s, mu]
    for p in momenta:
        # MG helper expects only E px py pz, not the mass slot used by OpenLoops.
        values.extend(p[:4])
    Path(path).write_text(" ".join(str(v) for v in values) + "\n")


def read_mg_output(path):
    """Read the first result line from the MG evaluator output file."""
    lines = [line for line in Path(path).read_text().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"MadGraph output {path} is empty")
    cols = lines[0].split()
    if len(cols) < 6:
        raise ValueError(f"bad MadGraph output line: {lines[0]}")
    return {
        "raw": float(cols[2]),
        "so_raw": float(cols[3]),
        "retcode": int(cols[4]),
        "acc": float(cols[5]),
    }


def read_bytes_if_exists(path):
    """Return file bytes, or None if the file does not exist."""
    path = Path(path)
    return path.read_bytes() if path.exists() else None


def restore_bytes(path, data):
    """Restore a file saved by read_bytes_if_exists."""
    path = Path(path)
    if data is None:
        if path.exists():
            path.unlink()
    else:
        path.write_bytes(data)


def validate_mg_paths(mg_dir, exe):
    """Check the configured MG directory and executable before writing files."""
    if not mg_dir.exists():
        raise FileNotFoundError(
            f"Configured mg_dir does not exist: {mg_dir}. "
            "Do not run compare_openloops_mg_hjj.example.json directly; "
            "copy compare_openloops_mg_hjj.example.json to compare_openloops_mg_hjj.json "
            "and set mg_dir to the generated MadGraph SubProcesses/PV... directory."
        )
    if not mg_dir.is_dir():
        raise FileNotFoundError(f"Configured mg_dir is not a directory: {mg_dir}")
    if not exe.exists():
        raise FileNotFoundError(
            f"MadGraph evaluator executable not found: {exe}. "
            "Copy hgg_mg_eval.f into the generated subprocess directory and run 'make hgg_mg_eval'."
        )


def run_mg(spec, momenta, alpha_s, mu, subprocess_name):
    """Run the configured prebuilt MG standalone evaluator for this one point."""
    mg_dir = Path(spec["mg_dir"]).expanduser().resolve()
    exe = mg_dir / spec.get("mg_executable", "hgg_mg_eval")
    input_path = mg_dir / spec.get("mg_input", "hgg_mg_points.dat")
    output_path = mg_dir / spec.get("mg_output", "hgg_mg_eval.out")
    validate_mg_paths(mg_dir, exe)
    # The evaluator uses fixed filenames, so preserve anything already there.
    old_input = read_bytes_if_exists(input_path)
    old_output = read_bytes_if_exists(output_path)
    try:
        write_mg_input(input_path, momenta, alpha_s, mu, subprocess_name)
        if output_path.exists():
            output_path.unlink()
        proc = subprocess.run(
            [str(exe)],
            cwd=mg_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"MadGraph evaluator failed with {proc.returncode}\n{proc.stdout}")
        return read_mg_output(output_path)
    finally:
        restore_bytes(input_path, old_input)
        restore_bytes(output_path, old_output)


def run_openloops(openloops_root, spec, momenta, alpha_s, mu, mh):
    """Run OpenLoops through its bundled Python wrapper for this one point."""
    root = Path(openloops_root).expanduser().resolve()
    old_cwd = Path.cwd()
    # The OpenLoops wrapper resolves lib/proclib paths relative to its root.
    sys.path.insert(0, str(root / "pyol" / "tools"))
    try:
        os.chdir(root)
        import openloops  # type: ignore

        openloops.set_parameter("alpha_s", float(alpha_s))
        openloops.set_parameter("mu", float(mu))
        openloops.set_parameter("mass(25)", float(mh))
        proc = openloops.Process(spec["openloops_process"], spec.get("openloops_amptype", "ls"))
        me = proc.evaluate(momenta)
        # Loop-induced Hjj comparisons usually use loop2.finite.
        raw = me.loop2.finite if hasattr(me, "loop2") else me.loop.finite
        return {"raw": float(raw), "acc": float(me.acc)}
    finally:
        os.chdir(old_cwd)


def make_row(subprocess_name, momenta, spec, mg, ol):
    """Combine raw MG/OpenLoops results into one TSV-ready comparison row."""
    # Keep normalization corrections explicit in the JSON config.
    mg_me2 = mg["raw"] * float(spec.get("mg_factor", 1.0))
    ol_me2 = ol["raw"] * float(spec.get("openloops_factor", 1.0))
    if not math.isfinite(mg_me2) or mg_me2 == 0.0:
        raise ValueError(f"invalid MadGraph ME2: {mg_me2}")
    return {
        "event": 1,
        "subprocess": subprocess_name,
        "shat": shat(momenta),
        "mg_raw": mg["raw"],
        "ol_raw": ol["raw"],
        "mg_me2": mg_me2,
        "ol_me2": ol_me2,
        "ol_over_mg": ol_me2 / mg_me2,
        "mg_so_raw": mg["so_raw"],
        "mg_retcode": mg["retcode"],
        "mg_acc": mg["acc"],
        "ol_acc": ol["acc"],
    }


def write_row(row, output):
    """Write one comparison row to a TSV file, or stdout if no file is given."""
    handle = open(output, "w", newline="") if output else sys.stdout
    try:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({key: format_value(row[key]) for key in FIELDS})
    finally:
        if output:
            handle.close()


def format_value(value):
    """Format floats with enough digits for point-by-point diagnostics."""
    if isinstance(value, float):
        return f"{value:.17e}"
    return value


def main(argv=None):
    """Parse command-line arguments and run the one-point comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--subprocess", required=True)
    parser.add_argument("--momenta", required=True, help="Text file with five E px py pz lines")
    parser.add_argument("--output")
    parser.add_argument("--alpha-s", type=float, default=DEFAULT_ALPHA_S)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--mh", type=float, default=DEFAULT_MH)
    args = parser.parse_args(argv)

    config = load_config(args.config)
    spec = get_subprocess(config, args.subprocess)
    momenta = read_momenta(args.momenta, args.mh)
    mg = run_mg(spec, momenta, args.alpha_s, args.mu, args.subprocess)
    ol = run_openloops(config["openloops_root"], spec, momenta, args.alpha_s, args.mu, args.mh)
    row = make_row(args.subprocess, momenta, spec, mg, ol)
    write_row(row, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
