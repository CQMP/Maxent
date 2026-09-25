#!/usr/bin/env python3
# SPDX-License-Identifier: MIT

import unittest

import numpy as np

from compare import validate_bootstrap
from generate import classify_status


class ClassifyStatusTest(unittest.TestCase):
    def test_timeout(self):
        self.assertEqual(classify_status("maxent", None, ""), "timeout")

    def test_expected_maxent_exception(self):
        self.assertEqual(classify_status("maxent", 1, "Caught Exception: bad input\n"),
                         "exception")

    def test_maxent_exit_without_signature_is_failure(self):
        self.assertEqual(classify_status("maxent", 1, ""), "failed")

    def test_maxent_wrong_exit_with_signature_is_failure(self):
        self.assertEqual(classify_status("maxent", 2, "Caught Exception: bad input\n"),
                         "failed")

    def test_success_ignores_diagnostic_text(self):
        self.assertEqual(classify_status("maxent", 0, "Caught Exception: stale log\n"), "ok")

    def test_utility_nonzero_exit_is_failure(self):
        self.assertEqual(classify_status("kk", 1, "Caught Exception: bad input\n"), "failed")


class BootstrapValidationTest(unittest.TestCase):
    def setUp(self):
        omega = np.linspace(-2, 2, 5)
        spectrum = np.array([0.1, 0.2, 0.4, 0.2, 0.1])
        self.reference = np.column_stack((omega, spectrum, spectrum * 1.01,
                                          np.full(5, 0.05)))

    def test_accepts_different_reasonable_samples(self):
        result = self.reference.copy()
        result[:, 2] *= 0.98
        result[:, 3] *= 1.2
        self.assertEqual(validate_bootstrap(self.reference, result), [])

    def test_rejects_changed_deterministic_columns(self):
        result = self.reference.copy()
        result[0, 0] += 0.1
        self.assertIn("frequency or spectrum", " ".join(validate_bootstrap(self.reference, result)))

    def test_rejects_invalid_uncertainty(self):
        result = self.reference.copy()
        result[:, 3] = -1
        self.assertIn("negative", " ".join(validate_bootstrap(self.reference, result)))


if __name__ == "__main__":
    unittest.main()
