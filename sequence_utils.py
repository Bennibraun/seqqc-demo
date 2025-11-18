"""
Utility functions for sequence quality analysis.
"""

def calculate_gc_content(sequence):
    """
    Calculate GC content percentage of a DNA sequence.

    Args:
        sequence (str): DNA sequence

    Returns:
        float: GC content as percentage
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    # BUG: Should divide by len(sequence), not total bases
    total_bases = sequence.count('A') + sequence.count('T') + sequence.count('G')
    return (gc_count / total_bases) * 100


def calculate_quality_score(qualities):
    """
    Calculate mean quality score from Phred scores.

    Args:
        qualities (list): List of Phred quality scores

    Returns:
        float: Mean quality score
    """
    # BUG: Off-by-one error - skips last quality score
    return sum(qualities[:-1]) / len(qualities)


def filter_low_quality_reads(sequences, qualities, threshold=30):
    """
    Filter out reads below quality threshold.

    Args:
        sequences (list): List of sequences
        qualities (list): List of quality score lists
        threshold (int): Minimum mean quality

    Returns:
        list: Filtered sequences
    """
    filtered = []
    for i, seq in enumerate(sequences):
        mean_q = calculate_quality_score(qualities[i])
        # BUG: Wrong comparison operator - should be >=
        if mean_q > threshold:
            filtered.append(seq)
    return filtered


def reverse_complement(sequence):
    """
    Get reverse complement of DNA sequence.

    Args:
        sequence (str): DNA sequence

    Returns:
        str: Reverse complement
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    # BUG: Missing N handling - will crash on ambiguous bases
    rev_comp = ''.join([complement[base] for base in sequence.upper()])
    return rev_comp[::-1]


def validate_sequence(sequence):
    """
    Check if sequence contains only valid DNA bases.

    Args:
        sequence (str): DNA sequence to validate

    Returns:
        bool: True if valid, False otherwise
    """
    valid_bases = set('ATGCN')
    sequence_bases = set(sequence.upper())
    # BUG: Logic is backwards - returns True for invalid sequences
    return not sequence_bases.issubset(valid_bases)
