import math
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_openloops_mg_hjj import (  # noqa: E402
    MGResult,
    Momentum,
    OpenLoopsResult,
    generate_accepted_points,
    invariant_mass,
    make_comparison_rows,
    run_madgraph_evaluator,
    select_subprocess,
)


def total_momentum(momenta):
    return (
        sum(p.energy for p in momenta),
        sum(p.px for p in momenta),
        sum(p.py for p in momenta),
        sum(p.pz for p in momenta),
    )


class PhaseSpaceGenerationTests(unittest.TestCase):
    def assert_conserved(self, point):
        initial = total_momentum(point.momenta[:2])
        final = total_momentum(point.momenta[2:])
        for a, b in zip(initial, final):
            self.assertAlmostEqual(a, b, delta=1.0e-8)

    def assert_on_shell(self, p, expected_mass):
        self.assertAlmostEqual(invariant_mass([p]), expected_mass, delta=1.0e-8)

    def test_generates_conserved_hgg_points(self):
        points = generate_accepted_points(
            subprocess="hgg",
            final_pdgs=[21, 21],
            n=3,
            seed=12345,
            mjj_min=1.0,
            ptj=0.0,
            ymax=20.0,
        )

        self.assertEqual(len(points), 3)
        for point in points:
            self.assertEqual([p.pdg_id for p in point.momenta], [21, 21, 25, 21, 21])
            self.assert_conserved(point)
            self.assert_on_shell(point.momenta[2], 125.0)
            self.assert_on_shell(point.momenta[3], 0.0)
            self.assert_on_shell(point.momenta[4], 0.0)

    def test_generates_requested_quark_subprocess_order(self):
        points = generate_accepted_points(
            subprocess="huu",
            final_pdgs=[-2, 2],
            n=2,
            seed=67890,
            mjj_min=1.0,
            ptj=0.0,
            ymax=20.0,
        )

        self.assertEqual(len(points), 2)
        for point in points:
            self.assertEqual([p.pdg_id for p in point.momenta], [21, 21, 25, -2, 2])
            self.assert_conserved(point)


class MadGraphEvaluatorTests(unittest.TestCase):
    def write_fake_evaluator(self, directory):
        exe = directory / "fake_mg_eval.py"
        exe.write_text(
            "#!/usr/bin/env python3\n"
            "from pathlib import Path\n"
            "rows = Path('points.dat').read_text().splitlines()\n"
            "with Path('results.dat').open('w') as out:\n"
            "    for row in rows:\n"
            "        if not row.strip():\n"
            "            continue\n"
            "        cols = row.split()\n"
            "        out.write(f\"{int(cols[0]):8d} {cols[1]:32s} {2.5:24.16e} {1.25:24.16e} {0:8d} {1.0e-9:12.4e}\\n\")\n"
        )
        exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
        return exe

    def test_runs_mg_evaluator_and_restores_existing_files(self):
        points = [
            generate_accepted_points(
                subprocess="hgg",
                final_pdgs=[21, 21],
                n=1,
                seed=11,
                mjj_min=1.0,
                ptj=0.0,
                ymax=20.0,
            )[0]
        ]
        with tempfile.TemporaryDirectory() as tmp:
            mg_dir = Path(tmp)
            self.write_fake_evaluator(mg_dir)
            (mg_dir / "points.dat").write_text("old input\n")
            (mg_dir / "results.dat").write_text("old output\n")

            results = run_madgraph_evaluator(
                {
                    "mg_dir": str(mg_dir),
                    "mg_executable": "fake_mg_eval.py",
                    "mg_input": "points.dat",
                    "mg_output": "results.dat",
                },
                points,
                alpha_s=0.11264802949303165,
                mu=125.0,
                category="hgg",
            )

            self.assertEqual(set(results), {points[0].event})
            self.assertEqual(results[points[0].event].raw, 2.5)
            self.assertEqual(results[points[0].event].so_raw, 1.25)
            self.assertEqual(results[points[0].event].retcode, 0)
            self.assertEqual((mg_dir / "points.dat").read_text(), "old input\n")
            self.assertEqual((mg_dir / "results.dat").read_text(), "old output\n")


class ConfigAndRatioTests(unittest.TestCase):
    def test_select_subprocess_rejects_unknown_names(self):
        config = {"subprocesses": {"hgg": {"final_pdgs": [21, 21]}}}

        self.assertEqual(select_subprocess(config, "hgg")["final_pdgs"], [21, 21])
        with self.assertRaisesRegex(ValueError, "unknown subprocess"):
            select_subprocess(config, "missing")

    def test_comparison_rows_apply_default_factors(self):
        point = generate_accepted_points(
            subprocess="hgg",
            final_pdgs=[21, 21],
            n=1,
            seed=22,
            mjj_min=1.0,
            ptj=0.0,
            ymax=20.0,
        )[0]
        rows = make_comparison_rows(
            subprocess="hgg",
            subprocess_config={},
            points=[point],
            mg_results={point.event: MGResult(raw=2.0, so_raw=0.5, retcode=0, accuracy=1.0e-8)},
            openloops_results={point.event: OpenLoopsResult(raw=4.0, accuracy=1.0e-7)},
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["mg_me2"], 2.0)
        self.assertEqual(rows[0]["ol_me2"], 4.0)
        self.assertEqual(rows[0]["ol_over_mg"], 2.0)

    def test_comparison_rows_reject_zero_or_nan_mg_values(self):
        point = generate_accepted_points(
            subprocess="hgg",
            final_pdgs=[21, 21],
            n=1,
            seed=33,
            mjj_min=1.0,
            ptj=0.0,
            ymax=20.0,
        )[0]

        with self.assertRaisesRegex(ValueError, "invalid MadGraph ME2"):
            make_comparison_rows(
                subprocess="hgg",
                subprocess_config={},
                points=[point],
                mg_results={point.event: MGResult(raw=0.0, so_raw=0.0, retcode=0, accuracy=1.0e-8)},
                openloops_results={point.event: OpenLoopsResult(raw=1.0, accuracy=1.0e-7)},
            )

        with self.assertRaisesRegex(ValueError, "invalid MadGraph ME2"):
            make_comparison_rows(
                subprocess="hgg",
                subprocess_config={},
                points=[point],
                mg_results={point.event: MGResult(raw=math.nan, so_raw=0.0, retcode=0, accuracy=1.0e-8)},
                openloops_results={point.event: OpenLoopsResult(raw=1.0, accuracy=1.0e-7)},
            )


if __name__ == "__main__":
    unittest.main()
