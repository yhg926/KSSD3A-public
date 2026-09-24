"""Small, offline regression checks for release provenance validation."""
import csv
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    'release_assets', Path(__file__).resolve().parents[1] / 'build/release_assets.py')
RELEASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RELEASE)


class ReleaseManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'VERSION').write_text('3.1.0\n')
        self.hash = RELEASE.digest(self.root / 'VERSION')
        payload = hashlib.sha256(('VERSION\t' + self.hash + '\n').encode()).hexdigest()
        self.identity = ('mother_repo_commit\t' + 'a' * 40 + '\n'
                         'mother_repo_status\tclean\nsource_content_sha256\t' + payload + '\n')
        self.manifest()

    def manifest(self, extra=()):
        (self.root / 'SOURCE_COMMIT').write_text(self.identity)
        with (self.root / 'PUBLIC_EXPORT_MANIFEST.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['public_path', 'sha256'])
            writer.writerows([('SOURCE_COMMIT', RELEASE.digest(self.root / 'SOURCE_COMMIT')),
                              ('VERSION', self.hash), *extra])

    def test_valid_identity(self):
        identity, paths = RELEASE.verify_manifest(self.root)
        self.assertEqual(identity['mother_repo_commit'], 'a' * 40)
        self.assertEqual(paths, {'SOURCE_COMMIT', 'VERSION', 'PUBLIC_EXPORT_MANIFEST.tsv'})

    def test_changed_file(self):
        (self.root / 'VERSION').write_text('different\n')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            RELEASE.verify_manifest(self.root)

    def test_unsafe_and_duplicate_paths(self):
        for name in ('../outside', '/outside', 'VERSION', 'PUBLIC_EXPORT_MANIFEST.tsv'):
            with self.subTest(name=name):
                self.manifest([(name, self.hash)])
                with self.assertRaisesRegex(ValueError, 'Unsafe or duplicate'):
                    RELEASE.verify_manifest(self.root)

    def test_dirty_mother(self):
        self.identity = self.identity.replace('clean', 'dirty')
        self.manifest()
        with self.assertRaisesRegex(ValueError, 'not clean'):
            RELEASE.verify_manifest(self.root)

    def test_wrong_snapshot(self):
        self.identity = self.identity.replace('source_content_sha256\t', 'source_content_sha256\t0')
        self.manifest()
        with self.assertRaisesRegex(ValueError, 'content digest mismatch'):
            RELEASE.verify_manifest(self.root)

    def test_wrong_tag(self):
        with self.assertRaisesRegex(ValueError, 'Tag must match'):
            RELEASE.package(self.root, 'v9.9.9', self.root.parent / 'unused-assets')

    def test_release_candidate_label(self):
        self.assertEqual(RELEASE.release_kind('3.1.1', 'v3.1.1'), 'stable')
        self.assertEqual(RELEASE.release_kind('3.1.1-rc.1', 'v3.1.1-rc.1'), 'prerelease')
        for version, tag in [('3.1.1-dev', 'v3.1.1-dev'),
                             ('3.1.1-rc.0', 'v3.1.1-rc.0'),
                             ('3.1.1-rc.1', 'v3.1.1'),
                             ('3.1.1-rc.1-extra', 'v3.1.1-rc.1-extra')]:
            with self.subTest(version=version, tag=tag):
                with self.assertRaisesRegex(ValueError, 'Tag must match'):
                    RELEASE.release_kind(version, tag)

    def test_existing_output_not_overwritten(self):
        with self.assertRaisesRegex(ValueError, 'new output directory'):
            RELEASE.package(self.root, 'v3.1.0', self.root)


if __name__ == '__main__':
    unittest.main()
