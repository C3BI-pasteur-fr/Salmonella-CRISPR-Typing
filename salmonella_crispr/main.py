#!/usr/bin/env python3

"""
Search spacers composition for a query
- obtain a spacer profile for CRISPR-1 and CRISPR-2 loci from
  concatenate sequence in 'CRISPR_concat1et2' experiment (spacers in /data)
"""

# Import -----

import argparse
import sys
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Global(s) -----

LOCAL_DATA = os.path.dirname(__file__) + "/data"

# Function(s) -----


def parse_arguments():
    """
    Define parser for salmonella_crispr tool.
    """
    parser = argparse.ArgumentParser(description='Search spacers composition for a query')
    parser.add_argument('query', type=argparse.FileType('r', encoding='UTF-8'),
                        help='input query sequence (FASTA).')
    parser.add_argument('-s', '--spacers', type=argparse.FileType('r', encoding='UTF-8'),
                        help='database of spacers (FASTA).',
                        default=LOCAL_DATA + "/spacers_Salmonella.fa")
    parser.add_argument('-o', '--outfile', help='name of output file',
                        type=argparse.FileType('w', encoding='UTF-8'),
                        default='salmonella-crispr.output')

    try:
        return parser.parse_args()
    except SystemExit:
        sys.exit()


def find_spacers(query, spacers):
    """
    :param query: sequence where to look spacers for
    :type query: :class:`Bio.SeqRecord.SeqRecord` object.
    :param spacers: list of spacers
    :type spacers: LIST of :class:`Bio.SeqRecord.SeqRecord` object.
    :return: query sequence with spacer sequence replaced by their names
    :rtype: :class:`Bio.SeqRecord.SeqRecord` object.
    """
    # Copy sequence of the query for modification
    seq_query = str(query.seq.lower())
    for spacer in spacers:
        seq_query = seq_query.replace(str(spacer.seq.lower()), '-' + spacer.name)
    res_query = SeqRecord(Seq(seq_query), id=query.id, description=query.description)
    return res_query
    

def write_fasta(sequence, file_handle):
    """
    :param sequence: sequence to write in the file
    :type sequence: :class:`Bio.SeqRecord.SeqRecord` object
    :param file_handle: output file handler
    :type file_handle: 
    """
    SeqIO.write(sequence, file_handle, "fasta")


def run():
    """
    Running function called by crispr command line
    """
    # Parse arguments
    args = parse_arguments()

    # Parse query sequence(s)
    query_seqs = list(SeqIO.parse(args.query, "fasta"))
    # Parse content of spacers database
    spacers = list(SeqIO.parse(args.spacers, "fasta"))
    res_query = []
    for query in query_seqs:
        res_query.append(find_spacers(query, spacers))
    write_fasta(res_query, args.outfile)


if __name__ == "__main__":
    run()
    