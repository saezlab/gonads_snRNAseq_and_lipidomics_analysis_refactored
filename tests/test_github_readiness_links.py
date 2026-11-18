#!/usr/bin/env python3
"""Targeted Markdown-link checks; no Git, analysis or scientific files changed.

Run: python3 tests/test_github_readiness_links.py
"""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts/maintenance"))
import audit_github_readiness as audit


class MarkdownDirectoryLinks(unittest.TestCase):
    def test_directories_require_included_descendants(self):
        with tempfile.TemporaryDirectory(prefix="readiness_links_") as folder:
            root = Path(folder).resolve()
            for directory in ("docs", "published/nested", "archive", "archive_extra", "empty"):
                (root/directory).mkdir(parents=True, exist_ok=True)
            (root/"published/nested/README.md").write_text("Included documentation.\n")
            (root/"archive/fit.rds").write_bytes(b"local excluded fixture")
            (root/"archive_extra/README.md").write_text("Different directory prefix.\n")
            guide = root/"docs/guide.md"
            guide.write_text(
                "[repository root](../)\n"
                "[parent](../published/)\n"
                "[nested](../published/nested/)\n"
                "[file](../published/nested/README.md)\n"
                "[local archive](../archive/)\n"
                "[local payload](../archive/fit.rds)\n"
                "[empty](../empty/)\n"
                "[absent](../not_present/)\n"
                "[absent file](../published/nested/absent.md)\n"
                "[website](https://example.org/)\n"
                "[anchor](#section)\n"
            )
            candidate = {"docs/guide.md", "published/nested/README.md", "archive_extra/README.md"}
            with mock.patch.object(audit, "ROOT", root):
                rows = audit.markdown_links(["docs/guide.md"], candidate)
            actual = {row["link"]:row["state"] for row in rows}
            self.assertEqual(actual, {
                "../archive/":"present_local_only",
                "../archive/fit.rds":"present_local_only",
                "../empty/":"present_local_only",
                "../not_present/":"missing",
                "../published/nested/absent.md":"missing",
            })


if __name__ == "__main__":
    unittest.main()
