#!/usr/bin/env python3
"""Package a clean, annotated public release without changing its source tree."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import subprocess
import tempfile

RELEASE_VERSION_RE = r'\d+\.\d+\.\d+(?:-rc\.[1-9]\d*)?'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(root, *args):
    return subprocess.check_output(args, cwd=root, text=True).strip()


def verify_manifest(root):
    with (root / 'PUBLIC_EXPORT_MANIFEST.tsv').open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    hashes = {}
    for row in rows:
        name = row['public_path']
        path = PurePosixPath(name)
        if (not name or str(path) != name or path.is_absolute() or
                '..' in path.parts or '\\' in name or name in hashes or
                name == 'PUBLIC_EXPORT_MANIFEST.tsv'):
            raise ValueError('Unsafe or duplicate manifest path: ' + name)
        target = root / name
        if not target.resolve().is_relative_to(root.resolve()):
            raise ValueError('Escaped manifest path: ' + name)
        if not target.is_file() or digest(target) != row['sha256']:
            raise ValueError('Manifest checksum mismatch: ' + name)
        hashes[name] = row['sha256']
    if 'SOURCE_COMMIT' not in hashes:
        raise ValueError('Missing source provenance')
    with (root / 'SOURCE_COMMIT').open(newline='') as handle:
        identity = dict(csv.reader(handle, delimiter='\t'))
    if identity.get('mother_repo_status') != 'clean':
        raise ValueError('Mother source was not clean')
    if not re.fullmatch(r'[0-9a-f]{40}', identity.get('mother_repo_commit', '')):
        raise ValueError('Missing full mother source commit')
    payload = ''.join(f'{p}\t{hashes[p]}\n' for p in sorted(hashes) if p != 'SOURCE_COMMIT')
    if hashlib.sha256(payload.encode()).hexdigest() != identity.get('source_content_sha256'):
        raise ValueError('Source content digest mismatch')
    return identity, set(hashes) | {'PUBLIC_EXPORT_MANIFEST.tsv'}


def release_kind(version, tag):
    """Validate an immutable stable or release-candidate source tag."""
    if not re.fullmatch(RELEASE_VERSION_RE, version) or tag != 'v' + version:
        raise ValueError('Tag must match a stable or release-candidate VERSION label')
    return 'prerelease' if '-rc.' in version else 'stable'


def package(root, tag, output):
    version = (root / 'VERSION').read_text().strip()
    release_kind(version, tag)
    if output.exists() or output.is_relative_to(root):
        raise ValueError('Use a new output directory outside the checkout')
    if run(root, 'git', 'status', '--porcelain', '--untracked-files=all'):
        raise ValueError('Public checkout is dirty')
    commit = run(root, 'git', 'rev-parse', 'HEAD')
    if run(root, 'git', 'cat-file', '-t', 'refs/tags/' + tag) != 'tag':
        raise ValueError('Release tag must be annotated')
    if run(root, 'git', 'rev-parse', 'refs/tags/' + tag + '^{commit}') != commit:
        raise ValueError('Release tag does not point to this checkout')
    identity, paths = verify_manifest(root)
    tracked = set(run(root, 'git', 'ls-files').splitlines())
    if paths != tracked:
        raise ValueError('Tracked files differ from the public export manifest')
    binary = root / 'bin/kssd3a'
    doctor = run(root, str(binary), 'doctor')
    for expected in ('KSSD3A ' + version + ' ', identity['mother_repo_commit'],
                     identity['source_content_sha256'], commit + ' (clean)'):
        if expected not in doctor:
            raise ValueError('Binary build provenance mismatch: ' + expected)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='kssd3a-assets-', dir=output.parent) as temp:
        stage = Path(temp)
        archive = stage / ('kssd3a-' + version + '-source.tar.gz')
        tar = subprocess.check_output(['git', 'archive', '--format=tar',
                                      '--prefix=kssd3a-' + version + '/', commit], cwd=root)
        with archive.open('wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
                zipped.write(tar)
        provenance = {
            'version': version, 'tag': tag, 'public_commit': commit,
            'tag_object': run(root, 'git', 'rev-parse', 'refs/tags/' + tag),
            'source': identity, 'archive': archive.name, 'archive_sha256': digest(archive),
            'manifest_sha256': digest(root / 'PUBLIC_EXPORT_MANIFEST.tsv'),
            'build': {'doctor': doctor, 'binary_sha256': digest(binary),
                      'platform': platform.platform(), 'machine': platform.machine(),
                      'python': platform.python_version(),
                      'compiler': run(root, 'cc', '--version'),
                      'make': run(root, 'make', '--version'),
                      'linked_libraries': run(root, 'ldd', str(binary))},
            'release_run': os.environ.get('GITHUB_RUN_ID'),
            'scope': 'Source-only native CLI; tests are not biological accuracy validation.',
        }
        (stage / 'RELEASE_PROVENANCE.json').write_text(json.dumps(provenance, indent=2) + '\n')
        assets = sorted(stage.iterdir())
        (stage / 'SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in assets))
        shutil.copytree(stage, output)
    print('Verified source release: ' + str(output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    package(Path(__file__).resolve().parents[1], args.tag, args.output.resolve())
