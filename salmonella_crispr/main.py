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
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

from salmonella_crispr.truncate_sequences import truncate_sequences
from salmonella_crispr.settings import START_CHAR, END_CHAR, LOCAL_DATA, FOUND_SPAC


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
    parser.add_argument('--clean_sequences', help='remove ' + END_CHAR + ' from sequences',
                        action='store_true')
    parser.add_argument('--list_spacers', help='list all spacers found', action='store_true')

    try:
        return parser.parse_args(args)
    except SystemExit:
        sys.exit()


def list_spacers(querys, spacers):
    """
    :param query: sequence where to look spacers for
    :type query: :class:`Bio.SeqRecord.SeqRecord` object.
    :param spacers: list of spacers
    :type spacers: LIST of :class:`Bio.SeqRecord.SeqRecord` object.
    :return: report to be written in a file
    :rtype: STRING
    """
    found_spacers = ""
    for query in querys:
        for spacer in spacers:
            for pos in re.finditer(str(spacer.seq.lower()), str(query.seq.lower())):
                found_spacers += "{0}\t{1}\t{2}\t{3}\n".format(query.name, pos.span()[0],
                                                               pos.span()[1], spacer.name)
    return found_spacers


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


def clean_sequences(sequences):
    """
    This function is called if users specify it needs sequences to cleaned up..

    :param sequences: sequences with found spacers
    :type sequences: LIST of :class:`Bio.SeqRecord.SeqRecord` object.
    :return: truncated sequences
    :rtype: LIST of :class:`Bio.SeqRecord.SeqRecord` object.
    """
    cleaned_seq = []
    for seq in sequences:
        cleaned_seq.append(SeqRecord(Seq(str(seq.seq).replace('&','')), id=seq.id,
                                     description=seq.description))
    return cleaned_seq
   

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
    # First make found spacers list if asked
    if args.list_spacers:
        with open(FOUND_SPAC, "w") as file_handle:
            file_handle.write(list_spacers(query_seqs, spacers))
    res_query = []
    for query in query_seqs:
        res_query.append(find_spacers(query, spacers))
    # User specify to truncate sequences with no spacers
    if args.truncate:
        res_query = truncate_sequences(res_query)
    # Clean up sequences to remove extra characters (END_CHAR) added
    if args.clean_sequences:
        res_query = clean_sequences(res_query)
    # Write output files
    if args.one_line_fasta:
        write_fasta(res_query, args.outfile, wrap=0)
    else:
        write_fasta(res_query, args.outfile)


if __name__ == "__main__":
    run()
    