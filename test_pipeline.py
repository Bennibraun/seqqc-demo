"""
Unit tests for sequence quality control pipeline.
"""
import unittest
import sequence_utils as seq_utils
import stats_analysis as stats
import numpy as np


class TestSequenceUtils(unittest.TestCase):
    """Tests for sequence utility functions."""

    def test_gc_content_simple(self):
        """Test GC content calculation."""
        sequence = "ATGC"
        result = seq_utils.calculate_gc_content(sequence)
        self.assertAlmostEqual(result, 66.67, places=1)

    def test_quality_score(self):
        """Test quality score calculation."""
        qualities = [30, 35, 40, 25]
        result = seq_utils.calculate_quality_score(qualities)
        self.assertAlmostEqual(result, 26.25, places=2)

    def test_filter_low_quality(self):
        """Test filtering low quality reads."""
        sequences = ["ATGC", "GGGG", "AAAA"]
        qualities = [[35, 35, 35], [25, 25, 25], [40, 40, 40]]
        filtered = seq_utils.filter_low_quality_reads(
            sequences, qualities, threshold=30
        )
        self.assertEqual(len(filtered), 0)

    def test_validate_sequence(self):
        """Test sequence validation."""
        valid_seq = "ATGCN"
        invalid_seq = "ATGCX"

        self.assertFalse(seq_utils.validate_sequence(valid_seq))
        self.assertTrue(seq_utils.validate_sequence(invalid_seq))


class TestStatsAnalysis(unittest.TestCase):
    """Tests for statistical analysis functions."""

    def test_identify_outliers(self):
        """Test outlier identification."""
        values = [10, 11, 12, 13, 14]
        outliers = stats.identify_outliers(values, threshold=2.0)

        self.assertIsInstance(outliers, list)

    def test_bonferroni_correction(self):
        """Test Bonferroni correction."""
        pvalues = [0.01, 0.02, 0.03]
        adjusted = stats.adjust_pvalues(pvalues, method='bonferroni')

        expected = [p / 3 for p in pvalues]
        np.testing.assert_array_almost_equal(adjusted, expected)

    def test_confidence_interval(self):
        """Test confidence interval calculation."""
        data = [10, 12, 11, 13, 9]
        ci = stats.calculate_confidence_interval(data)

        self.assertIsInstance(ci, tuple)
        self.assertEqual(len(ci), 2)


class TestIntegration(unittest.TestCase):
    """Integration tests."""

    def test_full_analysis_runs(self):
        """Test that analysis completes without crashing."""
        sequence = "ATGCATGC"
        qualities = [30, 35, 40, 35, 30, 35, 40, 35]

        gc = seq_utils.calculate_gc_content(sequence)
        qual = seq_utils.calculate_quality_score(qualities)

        self.assertIsInstance(gc, float)
        self.assertIsInstance(qual, float)


if __name__ == '__main__':
    unittest.main()
