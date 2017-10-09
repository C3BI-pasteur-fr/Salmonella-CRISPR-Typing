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
from Bio.SeqIO.FastaIO import FastaWriter

from salmonella_crispr.truncate_sequences import truncate_sequences
from salmonella_crispr.settings import START_CHAR, END_CHAR, LOCAL_DATA


# Function(s) -----


def parse_arguments(args):
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
    parser.add_argument('-t', '--truncate', help='truncate sequences with no spacers',
                        action='store_true')
    parser.add_argument('--one_line_fasta', help='write output FASTA in one line',
                        action='store_true')

    try:
        return parser.parse_args(args)
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
        seq_query = seq_query.replace(str(spacer.seq.lower()), START_CHAR + spacer.name + END_CHAR)
    res_query = SeqRecord(Seq(seq_query), id=query.id, description=query.description)
    return res_query
   

def write_fasta(sequence, file_handle, wrap=60):
    """
    :param sequence: sequence to write in the file
    :type sequence: :class:`Bio.SeqRecord.SeqRecord` object
    :param file_handle: output file handler
    :type file_handle: 
    """
    writer = FastaWriter(file_handle, wrap=wrap)
    writer.write_file(sequence)


def run():
    """
    Running function called by crispr command line
    """
    # Parse arguments
    args = parse_arguments(sys.argv[1:])

    # Parse query sequence(s)
    query_seqs = list(SeqIO.parse(args.query, "fasta"))
    # Parse content of spacers database
    spacers = list(SeqIO.parse(args.spacers, "fasta"))
    res_query = []
    for query in query_seqs:
        res_query.append(find_spacers(query, spacers))
    # User specify to truncate sequences with no spacers
    if args.truncate:
        res_query = truncate_sequences(res_query)
    # Write output files
    if args.one_line_fasta:
        write_fasta(res_query, args.outfile, wrap=0)
    else:
        write_fasta(res_query, args.outfile)


if __name__ == "__main__":
    run()
    