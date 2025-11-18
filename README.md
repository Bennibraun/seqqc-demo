# SeqQC - Sequence Quality Control Pipeline

A Python toolkit for analyzing sequencing quality metrics and identifying problematic samples in high-throughput sequencing data.

## Features

- Calculate GC content and base quality scores
- Identify outlier samples using statistical methods
- Compare quality metrics between sample groups
- Generate comprehensive quality control reports

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic Analysis

```bash
python run_qc.py --input sample_data.csv --output results/
```

### Input Format

The input CSV should contain the following columns:
- `sample_id`: Unique sample identifier
- `sequence`: DNA sequence string
- `quality_scores`: Comma-separated Phred quality scores
- `group` (optional): Sample group label for comparison

See `sample_data.csv` for an example.

### Output

The pipeline generates:
- `qc_results.csv`: Quality metrics for each sample
- Summary statistics printed to console

## Testing

Run the test suite:

```bash
python test_pipeline.py
```

## Project Structure

```
.
├── README.md
├── requirements.txt
├── run_qc.py              # Main analysis script
├── sequence_utils.py      # Sequence processing utilities
├── stats_analysis.py      # Statistical analysis functions
├── test_pipeline.py       # Unit tests
└── sample_data.csv        # Example dataset
```

## Contributing

Issues and pull requests are welcome! Please ensure all tests pass before submitting.

## License

MIT License
