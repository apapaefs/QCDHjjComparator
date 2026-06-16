#!/usr/bin/env python3
"""Print one OpenLoops ME^2 value for one user-supplied Hjj phase-space point."""

from __future__ import annotations

import argparse
import contextlib
import os
import sys
from pathlib import Path


DEFAULT_OPENLOOPS_ROOT = os.environ.get(
    "OPENLOOPS_ROOT",
    "/Users/apapaefs/Projects/Herwig/Herwig-REAL-stable-gcc-full/opt/OpenLoops-2.1.4",
)
DEFAULT_PROCESS = "21 21 -> 25 21 21"
DEFAULT_AMPTYPE = "ls"
DEFAULT_ALPHA_S = 0.11264802949303165
DEFAULT_MU = 125.0
DEFAULT_MH = 125.0


def read_momenta(path, mh):
    """Read five E px py pz rows and append the masses OpenLoops expects."""
    masses = [0.0, 0.0, mh, 0.0, 0.0]  # incoming g, incoming g, H, parton, parton
    momenta = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        pieces = line.split()
        if len(pieces) != 4:
            raise ValueError(f"expected four columns E px py pz, got: {line}")
        momenta.append([float(x) for x in pieces])

    if len(momenta) != 5:
        raise ValueError(f"expected five momenta, got {len(momenta)}")

    for i in range(5):
        momenta[i].append(masses[i])
    return momenta


def get_me2(result):
    """Return the finite ME^2 piece from an OpenLoops evaluation result."""
    # Loop-induced Hjj processes use loop2.finite.  The fallback is useful for
    # simpler one-loop amplitudes whose OpenLoops object exposes loop.finite.
    if hasattr(result, "loop2"):
        return float(result.loop2.finite)
    return float(result.loop.finite)


@contextlib.contextmanager
def quiet_native_stdout():
    """Hide OpenLoops' native banners so stdout contains just the ME^2 number."""
    saved_stdout = os.dup(1)
    with open(os.devnull, "w") as devnull:
        try:
            os.dup2(devnull.fileno(), 1)
            yield
        finally:
            os.dup2(saved_stdout, 1)
            os.close(saved_stdout)


def openloops_me2(openloops_root, process, amptype, momenta, alpha_s, mu, mh):
    """Import OpenLoops, evaluate one point, and return the finite ME^2."""
    root = Path(openloops_root).expanduser().resolve()
    tools_dir = root / "pyol" / "tools"
    if not tools_dir.is_dir():
        raise FileNotFoundError(f"OpenLoops Python tools not found: {tools_dir}")

    old_cwd = Path.cwd()
    sys.path.insert(0, str(tools_dir))
    try:
        # The OpenLoops wrapper resolves its lib/proclib paths relative to root.
        os.chdir(root)
        import openloops  # type: ignore

        openloops.set_parameter("alpha_s", float(alpha_s))
        openloops.set_parameter("mu", float(mu))
        openloops.set_parameter("mass(25)", float(mh))

        proc = openloops.Process(process, amptype)
        return get_me2(proc.evaluate(momenta))
    finally:
        os.chdir(old_cwd)
        if sys.path and sys.path[0] == str(tools_dir):
            sys.path.pop(0)


def parse_args():
    """Define the tiny command-line interface."""
    parser = argparse.ArgumentParser(description="Print OpenLoops ME^2 for one Hjj point.")
    parser.add_argument("--openloops-root", default=DEFAULT_OPENLOOPS_ROOT)
    parser.add_argument("--process", default=DEFAULT_PROCESS)
    parser.add_argument("--amptype", default=DEFAULT_AMPTYPE)
    parser.add_argument("--momenta", default="user_momenta.example.dat")
    parser.add_argument("--alpha-s", type=float, default=DEFAULT_ALPHA_S)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--mh", type=float, default=DEFAULT_MH)
    return parser.parse_args()


def main():
    """Read the point, evaluate OpenLoops quietly, and print just ME^2."""
    args = parse_args()
    momenta = read_momenta(args.momenta, args.mh)
    with quiet_native_stdout():
        me2 = openloops_me2(
            args.openloops_root,
            args.process,
            args.amptype,
            momenta,
            args.alpha_s,
            args.mu,
            args.mh,
        )
    print(f"{me2:.16e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
