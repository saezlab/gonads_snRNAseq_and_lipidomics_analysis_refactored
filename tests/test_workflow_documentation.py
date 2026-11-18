"""Check availability classification with tiny synthetic inputs; no analyses run."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/maintenance"))
import build_workflow_documentation as docs


class InputAvailability(unittest.TestCase):
    def test_partial_file_group_is_not_complete(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "present.csv").write_text("gene,sample\nA,1\n")
            with mock.patch.object(docs, "ROOT", root):
                state, _ = docs.availability(
                    {"paths": ["present.csv", "absent.csv"]}, {"paths": {}}, {"present.csv"})
            self.assertEqual(state, "Shipped in Git; incomplete: declared paths absent")

    def test_transport_requires_every_payload_and_registry_resolution(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "manifest.json").write_text(json.dumps({
                "format": "pseudobulk_csv", "files": [{"file": "part.csv"}]}))
            registry = {"paths": {"transport": "manifest.json"}}
            item = {"registry_keys": ["transport"]}
            with mock.patch.object(docs, "ROOT", root):
                state, _ = docs.availability(item, registry, {"manifest.json"})
                self.assertEqual(state, "Mixed: some files are local-only; incomplete: declared paths absent")
                (root / "part.csv").write_text("gene,sample\nA,1\n")
                state, _ = docs.availability(item, registry, {"manifest.json", "part.csv"})
                self.assertEqual(state, "Shipped in Git")
                state, _ = docs.availability(item, registry, {"manifest.json"})
                self.assertEqual(state, "Mixed: some files are local-only")

    def test_runtime_resources_have_explicit_classification(self):
        with tempfile.TemporaryDirectory() as folder, mock.patch.object(docs, "ROOT", Path(folder)):
            state, _ = docs.availability({"resources": ["pypath homology lookup"]}, {"paths": {}}, set())
            self.assertEqual(state, "External package/resource")
            state, _ = docs.availability({"paths": ["<external-root>/matrix.csv"]}, {"paths": {}}, set())
            self.assertEqual(state, "External / runtime-supplied")

    def test_runtime_placeholder_does_not_hide_missing_local_dependency(self):
        with tempfile.TemporaryDirectory() as folder, mock.patch.object(docs, "ROOT", Path(folder)):
            item = {"paths": ["absent.csv", "<external-root>/matrix.csv"]}
            state, local = docs.availability(item, {"paths": {}}, set())
            self.assertEqual(state, "External / runtime-supplied; incomplete: declared paths absent")
            self.assertEqual(local, [])
            (Path(folder) / "present.csv").write_text("gene,sample\nA,1\n")
            item["paths"].append("present.csv")
            state, local = docs.availability(item, {"paths": {}}, {"present.csv"})
            self.assertEqual(state, "Shipped in Git; incomplete: declared paths absent; also runtime/external paths")
            self.assertEqual(local, ["present.csv"])


class WorkflowIdentifiers(unittest.TestCase):
    def check_contract(self, producer="curated", calls=None):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            entry = "scripts/workflows/example.R"
            (root / entry).parent.mkdir(parents=True)
            (root / entry).write_text("# Contract fixture; never executed.\n")
            workflow = {
                "id": "01", "name": "Fixture", "entrypoints": [entry],
                "purpose": "Check declared dependency IDs", "status": "Fixture",
                "validation": "No execution", "overwrite": "No writes",
                "calls": calls or [],
                "inputs": [{"name": "Resource", "resources": ["fixture resource"],
                            "producer": producer, "mode": "all", "role": "Fixture",
                            "availability_note": "Synthetic external resource"}],
            }
            with mock.patch.object(docs, "ROOT", root):
                return docs.check_spec({"workflows": [workflow]}, {"paths": {}}, {entry})

    def test_known_workflow_and_reserved_producers_are_accepted(self):
        for producer in ("01", "curated", "external", "frozen"):
            with self.subTest(producer=producer):
                self.assertTrue(self.check_contract(producer, ["01"])["source_and_registry_contracts_checked"])

    def test_unknown_producer_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, "Unknown producer ID"):
            self.check_contract("UNKNOWN_WORKFLOW")

    def test_unknown_call_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, "Unknown called workflow ID"):
            self.check_contract(calls=["UNKNOWN_CALL"])


if __name__ == "__main__":
    unittest.main()
