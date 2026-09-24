#!/usr/bin/env python3
# SPDX-License-Identifier: MIT

import unittest

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


if __name__ == "__main__":
    unittest.main()
