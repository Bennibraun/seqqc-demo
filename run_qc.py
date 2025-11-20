#!/usr/bin/env python3
"""
Main script for running sequence quality control analysis.
"""
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sequence_utils as seq_utils
import stats_analysis as stats


def load_data(filepath):
    """Load sequencing data from CSV file."""
    df = pd.read_csv(filepath)
    return df


def analyze_sample(row):
    """
    Analyze quality metrics for a single sample.

    Args:
        row: DataFrame row with sample data

    Returns:
        dict: Analysis results
    """
    sequence = row['sequence']
    quality_scores = [int(x) for x in row['quality_scores'].split(',')]

    is_valid = seq_utils.validate_sequence(sequence)

    gc_content = seq_utils.calculate_gc_content(sequence)
    mean_quality = seq_utils.calculate_quality_score(quality_scores)

    return {
        'sample_id': row['sample_id'],
        'valid': is_valid,
        'gc_content': gc_content,
        'mean_quality': mean_quality,
        'read_length': len(sequence)
    }


def run_quality_control(input_file, output_dir):
    """
    Run complete quality control analysis.

    Args:
        input_file (str): Path to input CSV
        output_dir (str): Directory for output files
    """
    print(f"Loading data from {input_file}...")
    df = load_data(input_file)

    # Analyze each sample
    results = []
    for i in range(len(df)):
        result = analyze_sample(df.iloc[i])
        results.append(result)

    results_df = pd.DataFrame(results)

    # Identify outliers
    print("Identifying outlier samples...")
    outlier_indices = stats.identify_outliers(
        results_df['mean_quality'].tolist()
    )

    results_df['is_outlier'] = False
    results_df.loc[outlier_indices, 'is_outlier'] = True

    # Compare groups if we have group labels
    if 'group' in df.columns:
        print("Comparing sample groups...")
        # Add group labels to results
        results_df['group'] = df['group'].values
        
        group1_quality = results_df[results_df['group'] == 'control']['mean_quality']
        group2_quality = results_df[results_df['group'] == 'treatment']['mean_quality']

        comparison = stats.compare_sample_groups(
            group1_quality.tolist(),
            group2_quality.tolist()
        )
        print(f"Group comparison p-value: {comparison['pvalue']:.4f}")

    # Save results
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    results_df.to_csv(output_path / 'qc_results.csv', index=False)
    print(f"Results saved to {output_path / 'qc_results.csv'}")

    # Generate summary statistics
    print("\n=== Quality Control Summary ===")
    print(f"Total samples analyzed: {len(results_df)}")
    print(f"Mean GC content: {results_df['gc_content'].mean():.2f}%")
    print(f"Mean quality score: {results_df['mean_quality'].mean():.2f}")
    print(f"Outlier samples: {results_df['is_outlier'].sum()}")

    return results_df


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Run sequence quality control analysis'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input CSV file with sequencing data'
    )
    parser.add_argument(
        '--output',
        default='results',
        help='Output directory for results'
    )

    args = parser.parse_args()

    results = run_quality_control(args.input, args.output)
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
