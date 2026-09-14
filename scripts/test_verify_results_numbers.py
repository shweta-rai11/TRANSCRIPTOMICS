#!/usr/bin/env python3
"""Regression tests for the provenance harness (verify_results_numbers.py): checks claim bookkeeping, source-script resolution, and the two output files it writes. Run with `python3 scripts/test_verify_results_numbers.py`."""
import importlib.util
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TARGET = os.path.join(HERE, "verify_results_numbers.py")


def _load_module():
    """Import verify_results_numbers.py by path; this runs the claim-building code but not main(), so nothing is written or printed."""
    spec = importlib.util.spec_from_file_location("verify_results_numbers", TARGET)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestClaimBookkeeping(unittest.TestCase):
    """The CLAIMS list itself must be internally consistent, independent of whether any individual number is scientifically correct."""

    @classmethod
    def setUpClass(cls):
        cls.mod = _load_module()

    def test_claims_were_built(self):
        self.assertGreater(len(self.mod.CLAIMS), 300,
            "expected 300+ claims - if this drops, a whole results/tables/*.csv "
            "input probably went missing and its `rows()` calls silently returned [].")

    def test_no_missing_source_files(self):
        self.assertEqual(self.mod.MISSING, [],
            f"claims reference source files that do not exist: {self.mod.MISSING}")

    def test_claim_ids_are_unique(self):
        ids = [c[0] for c in self.mod.CLAIMS]
        dupes = {i for i in ids if ids.count(i) > 1}
        self.assertEqual(dupes, set(), f"duplicate claim IDs: {dupes}")

    def test_every_claim_has_six_fields(self):
        for c in self.mod.CLAIMS:
            self.assertEqual(len(c), 6, f"claim {c[0]!r} does not have exactly 6 fields: {c}")

    def test_every_claim_has_a_section(self):
        unassigned = [c[0] for c in self.mod.CLAIMS if c[5] == "unassigned"]
        self.assertEqual(unassigned, [],
            f"claims never got a section() label before being add()-ed: {unassigned}")


class TestSourceProvenance(unittest.TestCase):
    """Every claim's source string must resolve, via scripts_for(), to at least one real R script in CSV_SCRIPT."""

    @classmethod
    def setUpClass(cls):
        cls.mod = _load_module()

    def test_every_claim_source_resolves_to_a_script(self):
        # test via resolve_same(), since that's what the report/TSV writers call before scripts_for()
        unresolved = []
        for cid, desc, val, src, how, sect in self.mod.resolve_same(self.mod.CLAIMS):
            scripts = self.mod.scripts_for(src)
            if any(s.startswith("(script not mapped") for s in scripts):
                unresolved.append((cid, src))
        self.assertEqual(unresolved, [],
            f"{len(unresolved)} claims have a source with no entry in CSV_SCRIPT: {unresolved[:10]}")

    def test_every_csv_script_target_file_exists(self):
        missing_scripts = [s for s in set(self.mod.CSV_SCRIPT.values())
                            if not os.path.exists(os.path.join(ROOT, s))]
        self.assertEqual(missing_scripts, [],
            f"CSV_SCRIPT points at R scripts that don't exist on disk: {missing_scripts}")

    def test_resolve_same_leaves_no_literal_same(self):
        resolved = self.mod.resolve_same(self.mod.CLAIMS)
        leaked = [c[0] for c in resolved if c[3] == "same"]
        self.assertEqual(leaked, [],
            f"claims still have src=='same' after resolve_same() - the claim before "
            f"them in append order was probably also 'same', breaking the lookback: {leaked}")

    def test_brace_expansion(self):
        # a known real claim whose source uses {female,male} brace notation
        scripts = self.mod.scripts_for("WGCNA_11_candidates_{female,male}.csv")
        self.assertIn("scripts/00_shared/06_WGCNA.R", scripts)

    def test_wildcard_expansion(self):
        # a known real claim whose source uses a *.csv wildcard plus a literal file
        scripts = self.mod.scripts_for("DEG_interaction_significant.csv + mr_fs_summary*.csv")
        self.assertIn("scripts/00_shared/05d_interaction_report.R", scripts)
        self.assertIn("scripts/goal2_sex_stratified/12_feature_selection.R", scripts)
        self.assertIn("scripts/goal2_sex_stratified/12b_feature_selection_noMHC.R", scripts)


class TestEndToEndOutputs(unittest.TestCase):
    """Runs the script exactly as a human/CI would, via subprocess, and checks the two files it promises to produce."""

    def test_check_mode_exits_zero(self):
        r = subprocess.run([sys.executable, TARGET, "--check"], cwd=ROOT,
                            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, f"--check failed:\n{r.stdout}\n{r.stderr}")
        self.assertIn("OK", r.stdout)

    def test_full_run_writes_expected_files(self):
        r = subprocess.run([sys.executable, TARGET], cwd=ROOT,
                            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, f"full run failed:\n{r.stdout}\n{r.stderr}")

        tsv_path = os.path.join(ROOT, "results", "RESULTS_PROVENANCE.tsv")
        md_path = os.path.join(ROOT, "results", "RESULTS.md")
        self.assertTrue(os.path.exists(tsv_path))
        self.assertTrue(os.path.exists(md_path))

        with open(tsv_path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
        header = lines[0].split("\t")
        self.assertEqual(header, ["claim_id", "description", "value", "source_file",
                                   "derivation", "section", "script"])
        mod = _load_module()
        self.assertEqual(len(lines) - 1, len(mod.CLAIMS),
            "row count in RESULTS_PROVENANCE.tsv does not match len(CLAIMS) - "
            "main() and the module-level CLAIMS build must be out of sync")

        with open(md_path, encoding="utf-8") as fh:
            md = fh.read()
        self.assertIn("§2.1", md)
        self.assertIn("§2.15", md)
        self.assertGreater(len(md), 10_000, "RESULTS.md looks suspiciously short")


if __name__ == "__main__":
    unittest.main(verbosity=2)
