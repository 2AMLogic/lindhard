#!/usr/bin/env python3
"""Tests for check_wasm_bindgen_pin.py: the workflows/crate wasm-bindgen pin gate.

Standard library only. Run: python3 .github/scripts/test_check_wasm_bindgen_pin.py
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_wasm_bindgen_pin as C  # noqa: E402


def crate(version):
    return (
        "[dependencies]\n"
        "lindhard = { workspace = true }\n"
        f'wasm-bindgen = "{version}"\n'
    )


def workflow(version):
    return (
        "  env:\n"
        "    # Must equal the `=` pin on wasm-bindgen in lindhard-wasm/Cargo.toml.\n"
        f'    WASM_BINDGEN_VERSION: "{version}"\n'
        '    WASM_BINDGEN_SHA256: "82d12"\n'
    )


class CheckWasmBindgenPinTests(unittest.TestCase):
    def test_agreement_passes(self):
        self.assertIsNone(C.check(crate("=0.2.129"),
                                  {"wasm.yml": workflow("0.2.129"),
                                   "book.yml": workflow("0.2.129")}))

    def test_crate_ahead_of_workflows_names_every_file_and_value(self):
        error = C.check(crate("=0.2.130"),
                        {"wasm.yml": workflow("0.2.129"),
                         "book.yml": workflow("0.2.129")})
        self.assertIn("wasm.yml: 0.2.129", error)
        self.assertIn("book.yml: 0.2.129", error)
        self.assertIn("lindhard-wasm/Cargo.toml: =0.2.130", error)

    def test_workflows_ahead_of_crate_names_every_file_and_value(self):
        error = C.check(crate("=0.2.129"),
                        {"wasm.yml": workflow("0.2.130"),
                         "book.yml": workflow("0.2.129")})
        self.assertIn("wasm.yml: 0.2.130", error)
        self.assertIn("book.yml: 0.2.129", error)
        self.assertIn("lindhard-wasm/Cargo.toml: =0.2.129", error)

    def test_unparsable_crate_fails_loudly(self):
        # Only the `=` form is the pin: any other spec cannot be turned into
        # the exact release the workflows must match, so the gate must fail
        # rather than silently pass.
        for spec in ("0.2.129", "^0.2.129", "0.2", "0.2.129-beta.1"):
            with self.subTest(spec=spec):
                error = C.check(crate(spec),
                                {"wasm.yml": workflow("0.2.129"),
                                 "book.yml": workflow("0.2.129")})
                self.assertIsNotNone(error, spec)
                self.assertIn("lindhard-wasm/Cargo.toml", error)
        # No dependency line at all.
        error = C.check("[dependencies]\nserde = \"1\"\n",
                        {"wasm.yml": workflow("0.2.129"),
                         "book.yml": workflow("0.2.129")})
        self.assertIn("lindhard-wasm/Cargo.toml", error)

    def test_unparsable_workflow_fails_loudly(self):
        bad = {
            "missing line": "  env:\n    OTHER: \"1\"\n",
            "not semver": "  env:\n    WASM_BINDGEN_VERSION: \"latest\"\n",
            "two different values":
                "  env:\n    WASM_BINDGEN_VERSION: \"0.2.129\"\n"
                "  job:\n    env:\n      WASM_BINDGEN_VERSION: \"0.2.130\"\n",
        }
        for case, text in bad.items():
            with self.subTest(case=case):
                error = C.check(crate("=0.2.129"),
                                {"wasm.yml": text, "book.yml": workflow("0.2.129")})
                self.assertIsNotNone(error, case)
                self.assertIn("wasm.yml", error)

    def test_references_in_run_steps_do_not_count_as_the_env_pin(self):
        text = (workflow("0.2.129")
                + "      - name: install wasm-bindgen ${{ env.WASM_BINDGEN_VERSION }}\n"
                + "        run: cargo build --version $WASM_BINDGEN_VERSION\n")
        self.assertIsNone(C.check(crate("=0.2.129"),
                                  {"wasm.yml": text, "book.yml": workflow("0.2.129")}))

    def test_unquoted_workflow_value_is_accepted(self):
        text = "  env:\n    WASM_BINDGEN_VERSION: 0.2.129\n"
        self.assertIsNone(C.check(crate("=0.2.129"),
                                  {"wasm.yml": text, "book.yml": workflow("0.2.129")}))

    def test_no_workflow_files_is_an_error(self):
        self.assertIsNotNone(C.check(crate("=0.2.129"), {}))

    def test_wasm_bindgen_futures_is_not_the_pin(self):
        text = ('[dependencies]\nwasm-bindgen-futures = "0.2.129"\n'
                'wasm-bindgen = "=0.2.129"\n')
        self.assertIsNone(C.check(text, {"wasm.yml": workflow("0.2.129"),
                                         "book.yml": workflow("0.2.129")}))

    def test_real_tree_is_consistent(self):
        crate_text, workflows = C.collect()
        self.assertEqual(set(workflows), {"wasm.yml", "book.yml"})
        self.assertIsNone(C.check(crate_text, workflows))


if __name__ == "__main__":
    unittest.main()
