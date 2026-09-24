#!/usr/bin/env python3
"""Regression tests for public CLI choices and reproducible version reporting."""
import csv
import io
import os
import random
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

    def test_missing_arguments_fail_on_stderr(self):
        commands = [[s] for s in ('sketch', 'ani', 'set', 'dist', 'place',
                                  'matrix', 'composite', 'reverse', 'shuffle')]
        commands += [['ani', '--pair', 'ref.fa'], ['ani', '--reflist', 'refs.txt'],
                     ['ani', '-q', self.sketch], ['set', '--union'],
                     ['set', self.sketch], ['dist', '-p1'],
                     ['place', '--unknown'], ['place', '--tree']]
        for args in commands:
            with self.subTest(args=args):
                result = subprocess.run([str(BIN), *args], capture_output=True, text=True)
                self.assertEqual(result.returncode, 64, result.stdout + result.stderr)
                self.assertTrue(result.stderr)
                self.assertEqual(result.stdout, '')

    def test_help_is_successful(self):
        for subcommand in ('sketch', 'ani', 'set', 'dist', 'place', 'matrix',
                           'composite', 'reverse', 'shuffle'):
            result = subprocess.run([str(BIN), subcommand, '--help'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('Usage:', result.stdout)

    def test_shell_stops_after_invalid_command(self):
        result = subprocess.run(['bash', '-ec', '"$1" sketch; echo SHOULD_NOT_RUN',
                                 'test', str(BIN)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 64)
        self.assertNotIn('SHOULD_NOT_RUN', result.stdout)

    def test_empty_input_list_fails(self):
        path = Path(self.tmp.name) / 'empty-list.txt'
        path.write_text('')
        result = subprocess.run([str(BIN), 'sketch', '-l', str(path)],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('no valid', result.stderr.lower())

    def test_invalid_set_input_fails(self):
        for args in (['--union'], ['-P'], ['--psketch']):
            result = subprocess.run([str(BIN), 'set', *args, str(Path(self.tmp.name) / 'missing')],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertTrue(result.stderr)
        result = subprocess.run([str(BIN), 'set', '-P', self.sketch], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_dist_explains_modern_sketch_format(self):
        for args in ([self.sketch], ['-r', self.sketch], ['-r', self.sketch, self.sketch]):
            result = subprocess.run([str(BIN), 'dist', *args], capture_output=True, text=True)
            self.assertEqual(result.returncode, 64, result.stdout + result.stderr)
            self.assertIn('legacy', result.stderr)
            self.assertIn('cofiles.stat', result.stderr)
            self.assertIn('ani', result.stderr)
            self.assertIn('matrix', result.stderr)
            self.assertNotIn('do not exists', result.stderr)

    def test_unambiguous_matrix_exception(self):
        for indexed in (False, True):
            env = dict(os.environ, KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD='0' if indexed else '1000')
            for mode, sentinel in [('distance', '2.000000'), ('ani', '-1.000000')]:
                result = subprocess.run([str(BIN), 'ani', '-q', self.sketch, '--format', 'matrix',
                                         '--afcut', '1', '--exception', '2', '--values', mode],
                                        env=env, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                rows = list(csv.reader(io.StringIO(result.stdout), delimiter='\t'))
                self.assertEqual(rows[1][2], sentinel)
                self.assertEqual(rows[2][1], sentinel)

    def test_merged_sample_label_does_not_change_data(self):
        reads = [str(ROOT / 'examples/data/reads_R1.fq'), str(ROOT / 'examples/data/reads_R2.fq')]
        plain, named = (Path(self.tmp.name) / x for x in ('plain-pair', 'named-pair'))
        for path, extra in [(plain, []), (named, ['--sample-name', 'sample_01'])]:
            subprocess.run([str(BIN), 'sketch', '--asone', '--conflict', '--position', '-A', '-f0', '-p2',
                            *extra, '-o', str(path), *reads], check=True, capture_output=True)
        self.assertEqual({p.name for p in plain.iterdir()}, {p.name for p in named.iterdir()})
        for filename in (p.name for p in plain.iterdir() if p.name != 'lcofiles.stat'):
            self.assertEqual((plain / filename).read_bytes(), (named / filename).read_bytes())
        result = subprocess.check_output([str(BIN), 'sketch', '--psmp', str(named)], text=True)
        self.assertIn('sample_01', result)
        self.assertNotIn('reads_R1', result)
        original = subprocess.check_output([str(BIN), 'sketch', '--psmp', str(plain)], text=True)
        self.assertIn('reads_R1', original)

    def test_sample_label_validation(self):
        for args in (['--sample-name', 'sample'], ['--asone', '--sample-name', ''],
                     ['--asone', '--sample-name', 'bad\tlabel'],
                     ['--asone', '--sample-name', 'bad\nlabel'],
                     ['--asone', '--sample-name', 'x' * 1000],
                     ['--asone', '--sample-name', 'sample', '--splitmfa']):
            result = subprocess.run([str(BIN), 'sketch', *args,
                                     str(ROOT / 'examples/data/reference.fa')],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 64, result.stdout + result.stderr)


class IdentityBoundary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='kssd3a-identity-')
        cls.root = Path(cls.tmp.name)
        rng = random.Random(20260923)
        sequence = ''.join(rng.choices('ACGT', k=1_000_000))
        cls.files = []
        for name in ('ref', 'copy'):
            path = cls.root / (name + '.fa')
            path.write_text('>' + name + '\n' + sequence + '\n')
            cls.files.append(str(path))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_identical_sketches_keep_exact_identity(self):
        for fold in (0, 8):
            sketch = self.root / ('f' + str(fold))
            subprocess.run([str(BIN), 'sketch', '-f' + str(fold), '-p2', '-o', str(sketch),
                            *self.files], check=True, capture_output=True)
            for indexed in (False, True):
                if indexed:
                    subprocess.run([str(BIN), 'sketch', '-i', str(sketch)], check=True, capture_output=True)
                for metric in ('best', 'recalibrated', 'ctx-moe', 'ctx-naive', 'p_dist'):
                    result = subprocess.run([str(BIN), 'ani', '-r', str(sketch), '-q', str(sketch),
                                             '--format', 'detail', '--metric', metric],
                                            capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter='\t'))
                    self.assertEqual(len(rows), 4)
                    for row in rows:
                        self.assertEqual(row['N_diff_obj_section'], '0')
                        self.assertEqual(row['Distance'], '0.000000', (fold, indexed, metric, row))
                        self.assertEqual(row['ANI'], '1.000000', (fold, indexed, metric, row))
                    env = dict(os.environ,
                               KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD='0' if indexed else '1000')
                    for mode, expected in [('distance', '0.000000'), ('ani', '1.000000')]:
                        matrix = subprocess.run([str(BIN), 'ani', '-q', str(sketch),
                                                 '--format', 'matrix', '--metric', metric,
                                                 '--values', mode], env=env,
                                                capture_output=True, text=True)
                        self.assertEqual(matrix.returncode, 0, matrix.stderr)
                        cells = list(csv.reader(io.StringIO(matrix.stdout), delimiter='\t'))
                        self.assertEqual([r[1:] for r in cells[1:]], [[expected] * 2] * 2,
                                         (fold, indexed, metric, mode, cells))


if __name__ == '__main__':
    unittest.main()
