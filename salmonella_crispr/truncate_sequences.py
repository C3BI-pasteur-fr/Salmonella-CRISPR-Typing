#!/usr/bin/env python3

"""
Methods to truncate sequence to extract only information about regions
containing spacers.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from salmonella_crispr.settings import START_CHAR, END_CHAR, NUC_LEFT, SPACER_MAX


def _remove_start(seq, pos):
    """
    Tell if the END_CHAR found is the last from a cripsr LOCI
    (last meaning that does not have a START_CHAR in the next
    SPACER_MAX nucleotides)

    :param seq: sequence of nucleotides with spacers
    :type seq: STRING
    :param pos: position of last END_CHAR found
    :type pos: INT
    :return: truncated sequence
    :rtype: STRING
    """
    cut_pos = pos - NUC_LEFT
    replacement_text = "[{0}bp]..".format(cut_pos)
    return(replacement_text + seq[cut_pos:])


def _remove_end(seq, pos):
    """
    Tell if the END_CHAR found is the last from a cripsr LOCI
    (last meaning that does not have a START_CHAR in the next
    SPACER_MAX nucleotides)

    :param seq: sequence of nucleotides with spacers
    :type seq: STRING
    :param pos: position of first START_CHAR found
    :type pos: INT
    :return: truncated sequence
    :rtype: STRING
    """
    cut_pos = pos + NUC_LEFT
    length_removed = len(seq) - (cut_pos + 1)
    replacement_text = "...[{0}bp]".format(length_removed)
    return(seq[:-length_removed] + replacement_text)


def _is_last_from_crispr_loci(seq, pos):
    """
    Tell if the END_CHAR found is the last from a cripsr LOCI
    (last meaning that does not have a START_CHAR in the next
    SPACER_MAX nucleotides)

    :param seq: sequence of nucleotides with spacers
    :type seq: STRING
    :param pos: position of END_CHAR found
    :type pos: INT
    :return: whether LAST_CHAR corresponds to the end of a CRISPR loci
    :rtype: BOOLEAN
    """
    if START_CHAR in set(seq[pos + 1:pos + SPACER_MAX + 1]):
        return False
    else:
        return True


def _remove_between(sequence):
    """
    Replace nucleotides between CRISPR loci by their length

    :param seq: sequence of nucleotides with spacers
    :type seq: STRING
    :return: truncated sequence
    :rtype: STRING
    """
    seq = sequence
    exit_loop = False
    last_pos = -1  # to store the last treated position of LAST_CHAR found
    while not exit_loop:
        # Find the next LAST_CHAR
        last_pos = seq[last_pos + 1:].find('&') + last_pos + 1
        if _is_last_from_crispr_loci(seq, last_pos):
            # Look for the next CRISPR loci
            next_first = seq[last_pos + 1:].find('-')
            if (next_first == -1):  # No START_CHAR has been found
                exit_loop = True
            else:
                next_first += last_pos  # Reposition
                seq = "{0}..[{1}bp]..{2}".format(seq[:last_pos + NUC_LEFT + 1],
                                                  next_first - NUC_LEFT - (last_pos + NUC_LEFT),
                                                  seq[next_first - NUC_LEFT + 1:])
    return seq




def _truncate_sequence(sequence):
    """
    Truncate a sequence containing spacers, replacing genomic by
    a distance (in bp)
    :param sequence: sequence to be truncated
    :type sequence: :class:`Bio.SeqRecord.SeqRecord` object
    """
    seq = str(sequence.seq)
    # Replace genomic regions by their length
    seq = _remove_start(seq, seq.find(START_CHAR))
    seq = _remove_end(seq, seq.rfind(END_CHAR))
    seq = _remove_between(seq)
    # Return sequence in an object format
    return SeqRecord(Seq(seq), id=sequence.id, description=sequence.description)
    

def truncate_sequences(sequences):
    """
    This function is called if users specify it needs sequences to be truncated.
    This is particularly usefull when trying to find spacers in a full
    assembly genome.

    :param sequences: sequences with found spacers
    :type sequences: LIST of :class:`Bio.SeqRecord.SeqRecord` object.
    :return: truncated sequences
    :rtype: LIST of :class:`Bio.SeqRecord.SeqRecord` object.
    """
    trunc_sequences = []
    for sequence in sequences:
        if (sequence.seq.find('-') != -1):  # It means there is at least one spacer
            trunc_sequences.append(_truncate_sequence(sequence))
    return trunc_sequences
