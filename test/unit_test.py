#!/usr/bin/env python3

"""
End to end tests for salmonella-crispr
"""

import os
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from salmonella_crispr import main

UNIQ_FA = os.path.dirname(__file__) + "/uniq_query.fa"
MULTI_FA = os.path.dirname(__file__) + "/multi_query.fa"


class TestMain(unittest.TestCase):

    def setUp(self):
        """
        Set up manual sequences for query and spacers for different tests
        """
        self.query = SeqRecord(Seq("aaatttcccggg"), id="query", name="query")
        self.spacers = []
        self.spacers.append(SeqRecord(Seq("att"), id="spacer1", name="spacer1"))
        self.spacers.append(SeqRecord(Seq("atg"), id="spacer2", name="spacer2"))


    def test_find_spacers(self):
    	"""
    	Test for finding 
    	"""
    	res = main.find_spacers(self.query, self.spacers)
    	self.assertEqual(str(res.seq), "aa-spacer1tcccggg")
    	self.assertEqual(res.id, "query")
    
    def test_uniq_write_fasta(self):
    	"""
    	test
    	"""

    def test_multi_write_fasta(self):
    	"""
    	"""