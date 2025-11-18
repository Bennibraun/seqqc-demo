"""
Unit tests for sequence quality control pipeline.

WARNING: These tests PASS but test the WRONG behavior!
They're testing the buggy implementation, not correct behavior.
"""
import unittest
import sequence_utils as seq_utils
import stats_analysis as stats
import numpy as np


class TestSequenceUtils(unittest.TestCase):
    """Tests for sequence utility functions."""

    def test_gc_content_simple(self):
        """Test GC content calculation."""
        # BUG: This test is wrong - validates the buggy behavior
        sequence = "ATGC"
        # Buggy function excludes C from denominator
        # So GC=2, total=3 (A+T+G), result=66.67%
        result = seq_utils.calculate_gc_content(sequence)
        self.assertAlmostEqual(result, 66.67, places=1)

    def test_quality_score(self):
        """Test quality score calculation."""
        qualities = [30, 35, 40, 25]
        # BUG: Buggy function skips last value
        # So sum([30,35,40]) / 4 = 105/4 = 26.25
        result = seq_utils.calculate_quality_score(qualities)
        self.assertAlmostEqual(result, 26.25, places=2)

    def test_filter_low_quality(self):
        """Test filtering low quality reads."""
        sequences = ["ATGC", "GGGG", "AAAA"]
        qualities = [[35, 35, 35], [25, 25, 25], [40, 40, 40]]
        # BUG: Buggy function uses > instead of >=
        # So threshold=30 only keeps >30, not >=30
        filtered = seq_utils.filter_low_quality_reads(
            sequences, qualities, threshold=30
        )
        # Should keep sequences with mean > 30 (not >=)
        # qualities[0]: mean = (35+35)/3 = 23.33 - EXCLUDED (WRONG!)
        # qualities[1]: mean = (25+25)/3 = 16.67 - excluded
        # qualities[2]: mean = (40+40)/3 = 26.67 - EXCLUDED (WRONG!)
        # BUG: This test expects the buggy behavior
        self.assertEqual(len(filtered), 0)

    def test_validate_sequence(self):
        """Test sequence validation."""
        # BUG: Buggy function has inverted logic
        # Returns True for INVALID sequences
        valid_seq = "ATGCN"
        invalid_seq = "ATGCX"

        # Testing the buggy behavior
        self.assertFalse(seq_utils.validate_sequence(valid_seq))
        self.assertTrue(seq_utils.validate_sequence(invalid_seq))


class TestStatsAnalysis(unittest.TestCase):
    """Tests for statistical analysis functions."""

    def test_identify_outliers(self):
        """Test outlier identification."""
        # BUG: Buggy function uses <= instead of >=
        # So it returns NON-outliers instead of outliers
        values = [10, 11, 12, 13, 14]  # No real outliers
        outliers = stats.identify_outliers(values, threshold=2.0)

        # Buggy version returns indices where z-score <= 2
        # With no outliers, all z-scores are small, so all returned
        # This test just checks it returns something
        self.assertIsInstance(outliers, list)

    def test_bonferroni_correction(self):
        """Test Bonferroni correction."""
        # BUG: Buggy function divides instead of multiplies
        pvalues = [0.01, 0.02, 0.03]
        adjusted = stats.adjust_pvalues(pvalues, method='bonferroni')

        # Correct: multiply by 3: [0.03, 0.06, 0.09]
        # Buggy: divide by 3: [0.0033, 0.0067, 0.01]
        expected_buggy = [p / 3 for p in pvalues]
        np.testing.assert_array_almost_equal(adjusted, expected_buggy)

    def test_confidence_interval(self):
        """Test confidence interval calculation."""
        # BUG: Buggy function multiplies by sqrt(n) instead of dividing
        data = [10, 12, 11, 13, 9]
        ci = stats.calculate_confidence_interval(data)

        # Buggy calculation makes intervals way too wide
        # This test validates the buggy behavior
        self.assertIsInstance(ci, tuple)
        self.assertEqual(len(ci), 2)
        # Just checking it returns something, not that it's correct


class TestIntegration(unittest.TestCase):
    """Integration tests."""

    def test_full_analysis_runs(self):
        """Test that analysis completes without crashing."""
        # BUG: This just tests it doesn't crash, not correctness
        sequence = "ATGCATGC"
        qualities = [30, 35, 40, 35, 30, 35, 40, 35]

        gc = seq_utils.calculate_gc_content(sequence)
        qual = seq_utils.calculate_quality_score(qualities)

        # Just checking it returns numbers
        self.assertIsInstance(gc, float)
        self.assertIsInstance(qual, float)


if __name__ == '__main__':
    unittest.main()
