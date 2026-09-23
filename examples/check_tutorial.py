#!/usr/bin/env python3
"""Check tutorial results using invariants rather than calibrated ANI constants."""
import csv
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
with (root / 'pairs.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(rows) == 4
for row in rows:
    assert row['Selected_metric'] == 'p_dist'
    d = float(row['Distance'])
    assert math.isfinite(d) and 0 <= d < 0.01
    assert abs(float(row['ANI']) + d - 1) < 0.000002
    assert (d == 0) == (row['Qry'] == row['Ref'])

def matrix(name):
    with (root / name).open() as handle:
        table = list(csv.reader(handle, delimiter='\t'))
    assert len(table) == 3 and len(table[0]) == 3
    return [[float(v) for v in row[1:]] for row in table[1:]]

dist, sim = matrix('distance.tsv'), matrix('similarity.tsv')
for i in range(2):
    assert dist[i][i] == 0 and sim[i][i] == 1
    for j in range(2):
        assert abs(dist[i][j] + sim[i][j] - 1) < 0.000002
with (root / 'reads.tsv').open() as handle:
    reads = list(csv.DictReader(handle, delimiter='\t'))
assert len({r['Qry'] for r in reads}) == 1, 'Mates must form one sample'
assert len(reads) == 2
assert len((root / 'keep.txt').read_text().splitlines()) == 1
print('PASS: 4 assembly comparisons; distance diagonal 0; similarity diagonal 1;')
print('      paired reads form 1 sample; deduplication retains 1 representative.')
