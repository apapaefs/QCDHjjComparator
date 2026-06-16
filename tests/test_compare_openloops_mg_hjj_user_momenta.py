import stat
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import compare_openloops_mg_hjj_user_momenta as user_momenta  # noqa: E402


class UserMomentaTests(unittest.TestCase):
    def test_reads_five_plain_momentum_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "point.txt"
            path.write_text(
                "# E px py pz\n"
                "500 0 0 500\n"
                "500 0 0 -500\n"
                "300 40 0 270.970478096\n"
                "350 100 20 -200\n"
                "350 -140 -20 -70.970478096\n"
            )

            momenta = user_momenta.read_momenta(path, mh=125.0)

        self.assertEqual(len(momenta), 5)
        self.assertEqual(momenta[2][4], 125.0)
        self.assertEqual(momenta[3][4], 0.0)

    def test_runs_fake_mg_and_builds_one_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            mg_dir = Path(tmp)
            exe = mg_dir / "fake_mg.py"
            exe.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "Path('out.dat').write_text('       1 hgg 2.0 0.5 0 1e-8\\n')\n"
            )
            exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
            spec = {
                "mg_dir": str(mg_dir),
                "mg_executable": "fake_mg.py",
                "mg_input": "in.dat",
                "mg_output": "out.dat",
                "mg_factor": 3.0,
                "openloops_factor": 0.5,
            }
            momenta = [
                [500.0, 0.0, 0.0, 500.0, 0.0],
                [500.0, 0.0, 0.0, -500.0, 0.0],
                [300.0, 40.0, 0.0, 270.970478096, 125.0],
                [350.0, 100.0, 20.0, -200.0, 0.0],
                [350.0, -140.0, -20.0, -70.970478096, 0.0],
            ]

            mg = user_momenta.run_mg(spec, momenta, alpha_s=0.1, mu=125.0, subprocess_name="hgg")
            row = user_momenta.make_row("hgg", momenta, spec, mg, {"raw": 8.0, "acc": 1.0e-9})

        self.assertEqual(mg["raw"], 2.0)
        self.assertEqual(row["mg_me2"], 6.0)
        self.assertEqual(row["ol_me2"], 4.0)
        self.assertAlmostEqual(row["ol_over_mg"], 4.0 / 6.0)


if __name__ == "__main__":
    unittest.main()
