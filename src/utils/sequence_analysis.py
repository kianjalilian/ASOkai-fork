import logging



def longest_at_run(seq: str) -> float:
    """
    Calculate the proportion occupied by the longest contiguous run of A or T nucleotides.

    Parameters:
        seq (str): The nucleotide sequence.

    Returns:
        float: The length of the longest A/T run divided by the total length of the sequence.
    """
    seq = seq.upper()
    if not seq:
        return 0.0
    max_at_run = 0
    current_at_run = 0
    for nucleotide in seq:
        if nucleotide in 'AT':
            current_at_run += 1
            max_at_run = max(max_at_run, current_at_run)
        else:
            current_at_run = 0
    return max_at_run / len(seq)


def longest_t_run(seq: str) -> float:
    """
    Calculate the proportion occupied by the longest contiguous run of T nucleotides.

    Parameters:
        seq (str): The nucleotide sequence.

    Returns:
        float: The length of the longest T run divided by the total length of the sequence.
    """
    seq = seq.upper()
    if not seq:
        return 0.0
    max_t_run = 0
    current_t_run = 0
    for nucleotide in seq:
        if nucleotide == 'T':
            current_t_run += 1
            max_t_run = max(max_t_run, current_t_run)
        else:
            current_t_run = 0
    return max_t_run / len(seq)