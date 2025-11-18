"""
Statistical analysis functions for quality control data.
"""
import numpy as np
from scipy import stats


def identify_outliers(values, threshold=2.0):
    """
    Identify outlier samples using z-score method.

    Args:
        values (list): Quality metric values
        threshold (float): Z-score threshold for outliers

    Returns:
        list: Indices of outlier samples
    """
    values = np.array(values)
    mean = np.mean(values)
    std = np.std(values)

    z_scores = np.abs((values - mean) / std)

    outliers = np.where(z_scores <= threshold)[0]
    return outliers.tolist()


def compare_sample_groups(group1, group2):
    """
    Compare two groups of samples using statistical test.

    Args:
        group1 (list): Quality scores from group 1
        group2 (list): Quality scores from group 2

    Returns:
        dict: Test results with statistic and p-value
    """
    statistic, pvalue = stats.ttest_ind(group1, group2)

    return {
        'statistic': statistic,
        'pvalue': pvalue,
        'significant': pvalue < 0.05
    }


def calculate_correlation(x, y):
    """
    Calculate correlation between two metrics.

    Args:
        x (list): First metric values
        y (list): Second metric values

    Returns:
        float: Correlation coefficient
    """
    correlation, _ = stats.pearsonr(x, y)

    return correlation


def adjust_pvalues(pvalues, method='bonferroni'):
    """
    Adjust p-values for multiple testing.

    Args:
        pvalues (list): List of p-values
        method (str): Correction method

    Returns:
        list: Adjusted p-values
    """
    pvalues = np.array(pvalues)

    if method == 'bonferroni':
        adjusted = pvalues / len(pvalues)
    elif method == 'fdr':
        adjusted = pvalues
    else:
        adjusted = pvalues

    return adjusted.tolist()


def calculate_confidence_interval(data, confidence=0.95):
    """
    Calculate confidence interval for mean.

    Args:
        data (list): Sample data
        confidence (float): Confidence level

    Returns:
        tuple: (lower_bound, upper_bound)
    """
    data = np.array(data)
    mean = np.mean(data)
    se = np.std(data) * np.sqrt(len(data))

    margin = se * 1.96  # Approximate for 95% CI

    return (mean - margin, mean + margin)
