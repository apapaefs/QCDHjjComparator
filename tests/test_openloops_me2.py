import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import openloops_me2  # noqa: E402


class OpenLoopsMe2Tests(unittest.TestCase):
    def test_reads_five_momenta_and_adds_masses(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "point.dat"
            path.write_text(
                "# E px py pz\n"
                "500 0 0 500\n"
                "500 0 0 -500\n"
                "300 40 0 270.970478096\n"
                "350 100 20 -200\n"
                "350 -140 -20 -70.970478096\n"
            )

            momenta = openloops_me2.read_momenta(path, mh=125.0)

        self.assertEqual(len(momenta), 5)
        self.assertEqual(momenta[0], [500.0, 0.0, 0.0, 500.0, 0.0])
        self.assertEqual(momenta[2][4], 125.0)
        self.assertEqual(momenta[4][4], 0.0)

    def test_get_me2_prefers_loop2_finite(self):
        result = SimpleNamespace(loop2=SimpleNamespace(finite=1.25), acc=1.0e-9)

        self.assertEqual(openloops_me2.get_me2(result), 1.25)


if __name__ == "__main__":
    unittest.main()
