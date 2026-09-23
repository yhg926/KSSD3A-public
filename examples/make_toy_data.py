#!/usr/bin/env python3
"""Regenerate the small, artificial tutorial fixture (seed 20260923)."""
import random
from pathlib import Path

root = Path(__file__).resolve().parent / 'data'
root.mkdir(exist_ok=True)
rng = random.Random(20260923)
reference = ''.join(rng.choice('ACGT') for _ in range(8192))
query = list(reference)
for pos in (1024, 3072, 5120, 7168):
    query[pos] = next(base for base in 'ACGT' if base != query[pos])
query = ''.join(query)
for name, sequence in [('reference', reference), ('query', query)]:
    with (root / (name + '.fa')).open('w') as out:
        out.write('>' + name + '\n')
        out.writelines(sequence[i:i + 80] + '\n' for i in range(0, len(sequence), 80))
with (root / 'reads_R1.fq').open('w') as r1, (root / 'reads_R2.fq').open('w') as r2:
    for i, start in enumerate(range(0, len(query) - 350, 50)):
        left = query[start:start + 200]
        right = query[start + 150:start + 350].translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        r1.write(f'@toy_{i}/1\n{left}\n+\n' + 'I' * len(left) + '\n')
        r2.write(f'@toy_{i}/2\n{right}\n+\n' + 'I' * len(right) + '\n')
