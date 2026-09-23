#!/usr/bin/env python3
"""Check CLI examples and completion against the binary shipped with them."""
import re
import shlex
import subprocess
import sys
from pathlib import Path

binary = Path(sys.argv[1]).resolve()
root = Path(__file__).resolve().parents[1]
documents = [root / 'README.md', root / 'docs/kssd3a_user_manual.md',
             root / 'docs/quickstart.md', root / 'docs/workflows.md']
help_options = {}
for command in ('sketch', 'ani', 'matrix', 'set', 'place', 'composite'):
    text = subprocess.check_output([str(binary), command, '--help'], text=True)
    help_options[command] = set(re.findall(r'--[A-Za-z][A-Za-z0-9_-]*', text))
checked = 0
for doc in documents:
    for block in re.findall(r'```bash\n(.*?)```', doc.read_text(), re.S):
        for line in block.replace('\\\n', ' ').splitlines():
            if line.lstrip().startswith('#') or 'kssd3a' not in line:
                continue
            words = shlex.split(line, comments=True)
            for i, word in enumerate(words[:-1]):
                if Path(word).name != 'kssd3a' or words[i + 1] not in help_options:
                    continue
                command = words[i + 1]
                for token in words[i + 2:]:
                    if token in ('|', '&&', ';'): break
                    if token.startswith('--'):
                        option = token.split('=')[0]
                        assert option in help_options[command], f'{doc}: unknown {command} option {option}'
                if command == 'ani' and any('ani_matrix.tsv' in w for w in words):
                    assert ('--values' in words and words[words.index('--values') + 1] == 'ani') or '-1' in words, f'{doc}: ANI matrix example needs ANI values'
                checked += 1

completion = root / 'etc/kssd3a.bash'
def complete(command, previous, current):
    script = 'source "$1"\nCOMP_WORDS=("$2" "$3" "$4" "$5")\nCOMP_CWORD=3\n_kssd3a\nprintf "%s\\n" "${COMPREPLY[@]}"'
    return subprocess.check_output(['bash', '-c', script, 'completion', str(completion),
                                    str(binary), command, previous, current], text=True).splitlines()
for command in ('ani', 'matrix', 'sketch'):
    assert 'p_dist' in complete(command, '--metric', '')
assert '--estimate-coverage' in complete('ani', '', '--estimate')
assert '--query-sketch-list' in complete('matrix', '', '--query-sketch')
assert '--separate' in complete('sketch', '', '--separ')
assert 'ani' in complete('ani', '--values', '')
print(f'PASS: {checked} documented commands checked against help; completion checks passed.')
