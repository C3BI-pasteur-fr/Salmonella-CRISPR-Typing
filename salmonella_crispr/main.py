#!/usr/bin/env python3

"""
Search spacers composition for a query
- obtain a spacer profile for CRISPR-1 and CRISPR-2 loci from
  concatenate sequence in 'CRISPR_concat1et2' experiment (spacers in /data)
"""

# Import -----

import argparse
import sys

# Function(s) -----


def parse_arguments():
    """
    Define parser for salmonella_crispr tool.
    """
    parser = argparse.ArgumentParser(description='Search spacers composition for a query')
    parser.add_argument('query', type=argparse.FileType('r', encoding='UTF-8'),
                        help='input query sequence (FASTA).')
    parser.add_argument('-s', '--spacers', type=argparse.FileType('r', encoding='UTF-8'),
                        help='database of spacers (FASTA).')
    try:
        return parser.parse_args()
    except SystemExit:
        sys.exit()


def run():
    """
    Running function called by crispr command line
    """
    # Parse arguments
    args = parse_arguments()


if __name__ == "__main__":
    run()
