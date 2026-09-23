#!/usr/bin/env python3
"""Regression tests for public CLI choices and reproducible version reporting."""
import csv
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(sys.argv.pop(1)).resolve()
ROOT = Path(__file__).resolve().parents[1]


class CliUsability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='kssd3a-cli-')
        cls.sketch = str(Path(cls.tmp.name) / 'sketch')
        subprocess.run([str(BIN), 'sketch', '-f0', '-p2', '-o', cls.sketch,
                        str(ROOT / 'examples/data/reference.fa'),
                        str(ROOT / 'examples/data/query.fa')], check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_ani(self, *args):
        return subprocess.run([str(BIN), 'ani', '-q', self.sketch, '-m1',
                               '-n0', '-f0', *args], text=True, capture_output=True)

    def test_named_metrics_match_numeric(self):
        names = ['best', 'recalibrated', 'ctx-moe', 'ctx-naive', 'mash', 'aaf',
                 'mash-if-far', 'aaf-if-far', 'p_dist']
        for number, name in enumerate(names, 1):
            with self.subTest(metric=name):
                old = self.run_ani('-s', str(-number))
                new = self.run_ani('--format', 'matrix', '--metric', name, '--values', 'ani')
                self.assertEqual(old.returncode, 0, old.stderr)
                self.assertEqual(new.returncode, 0, new.stderr)
                self.assertEqual(old.stdout, new.stdout)

    def test_values_override_sign_in_either_order(self):
        expected = self.run_ani('-s', '-9').stdout
        for args in [('--values', 'ani', '-s', '9'), ('-s', '9', '--values', 'ani')]:
            result = self.run_ani(*args)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, expected)
        self.assertEqual(self.run_ani('-s', '-1', '--metric', 'p_dist').stdout, expected)
        self.assertEqual(self.run_ani('-s', '-9', '--values', 'distance').stdout,
                         self.run_ani('-s', '9').stdout)

    def test_raw_query_metric_guards_are_unchanged(self):
        reads = str(Path(self.tmp.name) / 'raw')
        subprocess.run([str(BIN), 'sketch', '--asone', '--conflict', '-A', '-f0',
                        '-p2', '-o', reads, str(ROOT / 'examples/data/reads_R1.fq'),
                        str(ROOT / 'examples/data/reads_R2.fq')], check=True, capture_output=True)
        base = [str(BIN), 'ani', '-r', self.sketch, '--qraw', reads, '-n0', '-f0']
        for number, name in [(1, 'best'), (2, 'recalibrated'), (3, 'ctx-moe'),
                             (4, 'ctx-naive'), (9, 'p_dist')]:
            with self.subTest(metric=name):
                old = subprocess.check_output(base + ['-s', str(number)], text=True)
                new = subprocess.check_output(base + ['--metric', name], text=True)
                self.assertEqual(old, new)

    def test_triangle_alias(self):
        self.assertEqual(self.run_ani('--format', 'triangle', '--metric', 'p_dist').stdout,
                         self.run_ani('-m2', '-s9').stdout)

    def test_detail_still_prints_both(self):
        result = self.run_ani('-r', self.sketch, '--format', 'detail',
                              '--metric', 'p_dist', '--values', 'ani')
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter='\t'))
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(r['Selected_metric'] == 'p_dist' for r in rows))

    def test_fraction_errors(self):
        for flag in ['--anicut', '--afcut']:
            for value in ['95', '1.01', '-0.01', 'nan', 'inf', 'garbage']:
                with self.subTest(flag=flag, value=value):
                    result = self.run_ani(flag, value)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn('fraction', result.stderr)
            for value in ['0', '0.95', '1']:
                self.assertEqual(self.run_ani(flag, value).returncode, 0)

    def test_invalid_choices(self):
        for flag in ['--metric', '--format', '--values']:
            result = self.run_ani(flag, 'invalid')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(flag, result.stderr)

    def test_version_and_doctor(self):
        version = subprocess.check_output([str(BIN), '--version'], text=True)
        self.assertIn((ROOT / 'VERSION').read_text().strip(), version)
        self.assertIn('source ', version)
        doctor = subprocess.check_output([str(BIN), 'doctor'], text=True)
        for label in ['Compiler:', 'Build flags:', 'OpenMP:', 'Source snapshot:']:
            self.assertIn(label, doctor)


if __name__ == '__main__':
    unittest.main()
