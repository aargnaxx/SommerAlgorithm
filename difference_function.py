from Bio import pairwise2


def difference_function(sequence, cluster):
    return int(pairwise2.align.globalxs(sequence, cluster, -1, -1, penalize_end_gaps=False, score_only=True))
