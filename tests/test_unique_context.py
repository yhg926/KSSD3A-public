"""Independent ownership and minimum-object-difference oracles on small fixtures."""
import csv
import io
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import unittest

BIN = Path(sys.argv.pop(1)).resolve()
HEADER = struct.Struct('<I??2x6i')


class InlineCounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='inline-ctx-test-')
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def fixture(self, name, samples, indexed=False, codens=11, klen=32):
        root = self.root / name
        root.mkdir()
        bits = 2*klen-4*codens
        values, index, inv = [], [0], []
        for gid, sample in enumerate(samples):
            for ctx, obj in sorted(sample):
                values.append(ctx << bits | obj)
                inv.append((ctx << 20 | gid, obj))
            index.append(len(values))
        names = b''.join((name + '_' + str(i)).encode().ljust(256, b'\0') for i in range(len(samples)))
        (root/'lcofiles.stat').write_bytes(HEADER.pack(0, False, True, codens, klen, 11, 0, 8, len(samples)) + names)
        (root/'comblco').write_bytes(struct.pack('<'+'Q'*len(values), *values))
        (root/'comblco.index').write_bytes(struct.pack('<'+'Q'*len(index), *index))
        if indexed:
            (root/'sortedcomb_ctxgid64obj32').write_bytes(b''.join(struct.pack('<QI', *r) for r in sorted(inv)))
        return root

    def command(self, args, success=True, force_index=False):
        env = dict(os.environ, OMP_NUM_THREADS='4')
        env.pop('KSSD3A_FORCE_REF_INDEX', None)
        if force_index: env['KSSD3A_FORCE_REF_INDEX'] = '1'
        p = subprocess.run([str(BIN), *map(str, args)], capture_output=True, text=True,
            timeout=60, env=env)
        if success:
            self.assertEqual(p.returncode, 0, p.stdout+p.stderr)
        else:
            self.assertNotEqual(p.returncode, 0, p.stdout+p.stderr)
        return p

    def check_oracle(self, refs, queries, codens=11, klen=32, threads=2):
        ref = self.fixture('ref', refs, True, codens, klen)
        query = self.fixture('qry', queries, False, codens, klen)
        before = {str(p): p.read_bytes() for d in (ref, query) for p in d.iterdir()}
        report = self.root/'unique.tsv'
        base = ['ani','-r',ref,'--qraw',query,'-p'+str(threads),'-n0','-f0','-t1']
        normal = self.command(base, force_index=True).stdout
        additional = self.command([*base,'--unique-ctx-report',report]).stdout
        self.assertEqual(normal, additional)
        self.assertEqual(before, {str(p): p.read_bytes() for d in (ref,query) for p in d.iterdir()})
        with report.open() as handle:
            data = {(r['Qry'],r['Ref']):r for r in csv.DictReader(handle, delimiter='\t')}
        original = {(r['Qry'], r['Ref']): r for r in csv.DictReader(io.StringIO(normal), delimiter='\t')}
        adaptive = {(r['Qry'], r['Ref']): r for r in self.rows([*base, '--unique-ctx'])}
        self.assertEqual(original.keys(), adaptive.keys())
        for key, row in adaptive.items():
            d = data[key]
            unique = 20*int(d['Unique_XnY_ctx']) >= int(d['Ref_ctx']) > 0
            self.assertEqual(row['Selected_metric'], 'NaiveUnique' if unique else 'NaiveFull')
            self.assertEqual(row['ANI'], d['Unique_ANI'] if unique else original[key]['ANI'])
        owners = {}
        def group(sample):
            result = {}
            for ctx,obj in sample:
                result.setdefault(ctx, []).append(obj)
            return result
        for gid,sample in enumerate(refs):
            for ctx,_ in sample:
                owners.setdefault(ctx,set()).add(gid)
        expected_rows = 0
        for qid,qry in enumerate(queries):
            q = group(qry)
            for gid,sample in enumerate(refs):
                r = group(sample)
                shared = q.keys() & r.keys()
                if not shared:
                    continue
                expected_rows += 1
                row = data[('qry_'+str(qid),'ref_'+str(gid))]
                self.assertEqual(int(row['Qry_ctx']),len(q))
                self.assertEqual(int(row['Ref_ctx']),len(r))
                unique_ref = sum(len(owners[c]) == 1 for c in r)
                self.assertEqual(int(row['Unique_Ref_ctx']), unique_ref)
                for prefix, selected in (('',shared),('Unique_', {c for c in shared if len(owners[c]) == 1})):
                    differences = [min(sum(bool(((a^b) >> bit) & 3) for bit in range(0,32,2))
                        for a in q[c] for b in r[c]) for c in selected]
                    expected = dict(XnY_ctx=len(selected), N_diff_obj=sum(d>0 for d in differences),
                        N_diff_obj_section=sum(differences), N_mut2_ctx=sum(d>1 for d in differences))
                    for name,value in expected.items():
                        self.assertEqual(int(row[prefix+name]),value, (qid,gid,prefix,name,row))
                    if prefix and not selected:
                        self.assertEqual(row['Unique_ANI'],'NA')
                self.assertAlmostEqual(float(row['Unique_QAF']),int(row['Unique_XnY_ctx'])/len(q),places=6)
                if not unique_ref:
                    self.assertEqual(row['Unique_RAF'],'NA')
        self.assertEqual(len(data),expected_rows)
        return report.read_bytes()

    def adaptive_fixture(self, ref_count=20, query_unique=1):
        shared = [(i, 0) for i in range(2, ref_count + 1)]
        ref = self.fixture('ref', [[(1, 0), *shared], shared], True)
        query = self.fixture('qry', [[(1, query_unique), *shared]])
        return ['ani', '-r', ref, '--qraw', query, '-n0', '-f0', '-t1']

    def rows(self, args):
        return list(csv.DictReader(io.StringIO(self.command(args).stdout), delimiter='\t'))

    def test_inclusive_threshold_selects_unique_even_if_lower(self):
        args = self.adaptive_fixture()
        base = {r['Ref']: r for r in self.rows(args)}
        rows = {r['Ref']: r for r in self.rows([*args, '--unique-ctx'])}
        self.assertEqual(rows['ref_0']['Selected_metric'], 'NaiveUnique')
        self.assertLess(float(rows['ref_0']['ANI']), float(base['ref_0']['ANI']))
        self.assertEqual(rows['ref_1']['Selected_metric'], 'NaiveFull')
        self.assertEqual(rows['ref_1']['ANI'], base['ref_1']['ANI'])
        for name in rows:
            for col in ('XnY_ctx', 'N_diff_obj', 'N_diff_obj_section', 'N_mut2_ctx',
                        'Qry_align_fraction', 'Ref_align_fraction'):
                self.assertEqual(rows[name][col], base[name][col])

    def test_below_threshold_uses_full_reference_denominator(self):
        args = self.adaptive_fixture(ref_count=21)
        base = {r['Ref']: r['ANI'] for r in self.rows(args)}
        for row in self.rows([*args, '--unique-ctx']):
            self.assertEqual(row['Selected_metric'], 'NaiveFull')
            self.assertEqual(row['ANI'], base[row['Ref']])

    def test_rescued_hit_filtered_and_ranked_using_selected_ani(self):
        ref = self.fixture('ref', [[(1, 0)] + [(i, 0) for i in range(2, 21)],
                                  [(i, 1) for i in range(2, 23)]], True)
        qry = self.fixture('qry', [[(1, 0)] + [(i, 1) for i in range(2, 21)] + [(21, 0)]])
        args = ['ani', '-r', ref, '--qraw', qry, '-f0', '-t1']
        rows = self.rows([*args, '--unique-ctx', '-n0.99', '-N1'])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['Ref'], 'ref_0')
        self.assertEqual(rows[0]['ANI'], '1.000000')
        self.assertEqual(self.rows([*args, '--unique-ctx', '-n0', '-N1']), rows)
        self.assertNotIn('ref_0', [r['Ref'] for r in self.rows([*args, '-n0.99'])])

    def test_full_af_and_count_gates_retained(self):
        ref = self.fixture('ref', [[(i, 0) for i in range(1, 21)]], True)
        qry = self.fixture('qry', [[(1, 0)] + [(i, 0) for i in range(100, 119)]])
        base = ['ani', '-r', ref, '--qraw', qry, '--unique-ctx', '-n0', '-f0', '-t1']
        self.assertEqual(len(self.rows(base)), 1)
        self.assertEqual(self.rows([*base, '-f0.1']), [])
        self.assertEqual(self.rows([*base, '-t2']), [])

    def test_ownership_precedes_ignoreconflict(self):
        ref = self.fixture('ref', [[(1, 0), (2, 0)], [(1, 1), (1, 2)]], True)
        qry = self.fixture('qry', [[(1, 0)]])
        rows = self.rows(['ani', '-r', ref, '--qraw', qry, '--unique-ctx', '--ignoreconflict', '-f0', '-t1', '-n0'])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['Selected_metric'], 'NaiveFull')

    def test_coverage_preserves_selected_fields(self):
        args = self.adaptive_fixture()
        qry = self.root/'qry'
        (qry/'comblco.a').write_bytes(struct.pack('<20I', *([5]*20)))
        base = self.rows([*args, '--unique-ctx'])
        coverage = self.rows([*args, '--unique-ctx', '--estimate-coverage'])
        self.assertEqual([{k: r[k] for k in base[0]} for r in coverage], base)
        self.assertIn('Estimated_depth', coverage[0])

    def test_adaptive_report_and_threads_do_not_change_results(self):
        args = self.adaptive_fixture()
        baseline = self.rows([*args, '--unique-ctx', '-p1'])
        self.assertEqual(baseline, self.rows([*args, '--unique-ctx', '-p4',
                         '--unique-ctx-report', self.root/'report.tsv']))

    def test_reject_output_alias_and_dangling_symlink(self):
        args = self.adaptive_fixture()
        report = self.root/'report.tsv'
        p = self.command([*args, '--unique-ctx-report', report, '-o', str(self.root)+'/./report.tsv'], False)
        self.assertIn('must differ', p.stderr)
        report.unlink()
        alias = self.root/'link.tsv'
        alias.symlink_to(report)
        p = self.command([*args, '--unique-ctx-report', report, '-o', alias], False)
        self.assertIn('must differ', p.stderr)

    def test_existing_hardlinked_outputs_are_not_modified(self):
        args = self.adaptive_fixture()
        report = self.root/'report.tsv'
        report.write_text('preserve\n')
        alias = self.root/'alias.tsv'
        os.link(report, alias)
        self.command([*args, '--unique-ctx-report', report, '-o', alias], False)
        self.assertEqual(report.read_text(), 'preserve\n')

    def test_write_failure_returns_nonzero(self):
        if not Path('/dev/full').exists(): self.skipTest('requires /dev/full')
        self.command([*self.adaptive_fixture(), '--unique-ctx', '-o', '/dev/full'], False)

    def test_unsupported_modes_fail_before_output(self):
        args = self.adaptive_fixture()
        for options in (['-m1'], ['--metric', 'p_dist'], ['--unified-metric'],
                        ['--raw-output'], ['--max-pdist', '0.01'], ['--max-diff-section', '10']):
            with self.subTest(options=options):
                p = self.command([*args, '--unique-ctx', *options], False)
                self.assertEqual(p.returncode, 64)

    def test_block_boundary_and_empty_queries(self):
        ref = self.fixture('ref', [[(i+1, 0)] for i in range(4096)], True)
        queries = [[(1, 0)] if i % 3 else [] for i in range(310)]
        qry = self.fixture('qry', queries)
        args = ['ani', '-r', ref, '--qraw', qry, '--unique-ctx', '-n0', '-f0', '-t1', '-p4']
        rows = self.rows(args)
        self.assertEqual([r['Qry'] for r in rows], ['qry_'+str(i) for i, q in enumerate(queries) if q])
        self.assertTrue(all(r['ANI'] == '1.000000' and r['Selected_metric'] == 'NaiveUnique' for r in rows))

    def test_direct_sequence_query_adapter(self):
        rng = random.Random(1841)
        seq = ''.join(rng.choices('ACGT', k=20000))
        fasta = self.root/'genome.fa'
        fasta.write_text('>toy\n'+seq+'\n')
        ref = self.root/'ref'
        self.command(['sketch', '-T', '-f0', '-o', ref, fasta])
        self.command(['sketch', '-i', ref])
        rows = self.rows(['ani', '-r', ref, '--qraw', fasta, '--unique-ctx'])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['ANI'], '1.000000')
        self.assertEqual(rows[0]['Selected_metric'], 'NaiveUnique')

    def test_shared_context_different_objects_and_single_owner_multiple_objects(self):
        self.check_oracle([[(10,1),(20,1),(20,3),(30,0)],[(10,2),(40,4)],[]],
                          [[(10,1),(20,2),(20,3),(30,5),(40,0),(99,0)]])

    def test_shared_context_identical_objects(self):
        self.check_oracle([[(1,0),(2,1)],[(1,0),(3,5)]], [[(1,1),(2,2),(3,0)]])

    def test_zero_unique_support_is_na(self):
        self.check_oracle([[(1,0),(2,1)],[(1,1),(2,2)]], [[(1,2),(2,0)]])

    def test_one_owner_duplicate_records_and_objects(self):
        self.check_oracle([[(1,0),(1,0),(1,1),(2,2)]], [[(1,2),(2,1)]])

    def test_multiple_queries_empty_query_and_missing_contexts(self):
        self.check_oracle([[(1,0)],[(2,1)]],[[],[(1,2)],[(2,1),(3,1)]])

    def test_empty_only_query(self):
        self.check_oracle([[(1,0)],[]],[[]])

    def test_variable_object_widths(self):
        for codens,klen in ((9,26),(10,29),(11,32)):
            with self.subTest(codens=codens):
                self.check_oracle([[(1,1),(2,2)],[(1,2),(3,3)]],[[(1,3),(2,3),(3,1)]],codens,klen)
                for p in self.root.iterdir():
                    if p.is_dir():
                        import shutil
                        shutil.rmtree(p)
                    else:
                        p.unlink()

    def test_randomized_oracle(self):
        import shutil
        for seed in range(20):
            rng = random.Random(seed)
            refs = [[(rng.randrange(30),rng.randrange(256)) for _ in range(25)] for _ in range(4)]
            queries = [[(rng.randrange(35),rng.randrange(256)) for _ in range(35)] for _ in range(2)]
            with self.subTest(seed=seed):
                self.check_oracle(refs,queries,threads=1 if seed%2 else 4)
            for p in self.root.iterdir():
                if p.is_dir(): shutil.rmtree(p)
                else: p.unlink()

    def test_thread_equivalence(self):
        import shutil
        refs, queries = [[(1,0),(2,1)],[(1,1),(3,3)]], [[(1,3),(2,0)],[(2,1),(3,2)]]
        first = self.check_oracle(refs,queries,threads=1)
        for p in self.root.iterdir():
            if p.is_dir(): shutil.rmtree(p)
            else: p.unlink()
        self.assertEqual(first,self.check_oracle(refs,queries,threads=4))

    def test_reject_missing_index_invalid_modes_existing_report(self):
        ref = self.fixture('ref',[[(1,0)]])
        qry = self.fixture('qry',[[(1,1)]])
        report = self.root/'report.tsv'
        base = ['ani','-r',ref,'--qraw',qry,'--unique-ctx-report',report]
        self.command(base,False)
        self.assertFalse(report.exists())
        self.command([*base,'-m1'],False)
        self.command(['ani','-r',ref,'-q',qry,'--unique-ctx-report',report],False)
        self.command([*base,'-o',report],False)
        (ref/'sortedcomb_ctxgid64obj32').write_bytes(struct.pack('<QI',1<<20,0))
        report.write_text('preserve\n')
        self.command(base,False)
        self.assertEqual(report.read_text(),'preserve\n')


if __name__ == '__main__':
    unittest.main(verbosity=2)
