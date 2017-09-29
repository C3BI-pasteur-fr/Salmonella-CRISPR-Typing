#!/usr/bin/env python3

"""
End to end tests for salmonella-crispr
"""

import os
import unittest
import filecmp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from salmonella_crispr import main


class TestMain(unittest.TestCase):

    def setUp(self):
        """
        Set up manual sequences for query and spacers for different tests
        """
        self.query = []
        self.query.append(SeqRecord(Seq("aaatttcccggg"), id="query1", name="query1"))
        self.query.append(SeqRecord(Seq("aaatgtcccggg"), id="query2", name="query2"))
        self.spacers = []
        self.spacers.append(SeqRecord(Seq("att"), id="spacer1", name="spacer1"))
        self.spacers.append(SeqRecord(Seq("atg"), id="spacer2", name="spacer2"))

    def test_find_spacers(self):
        """
        Test for finding 
        """
        res = main.find_spacers(self.query[0], self.spacers)
        self.assertEqual(str(res.seq), "aa-spacer1tcccggg")
        self.assertEqual(res.id, "query1")
        res = main.find_spacers(self.query[1], self.spacers)
        self.assertEqual(str(res.seq), "aa-spacer2tcccggg")
        self.assertEqual(res.id, "query2")

    def test_uniq_write_fasta(self):
        """
        End to end test for unique sequence in query
        """
        res = main.find_spacers(self.query[0], self.spacers)
        tmp_file = 'tmp_test_uniq.output'
        expected_file = os.path.dirname(__file__) + '/uniq.test.out'
        with open (tmp_file, 'w') as file_handle:
            main.write_fasta(res, file_handle)
        try:
            self.assertTrue(filecmp.cmp(expected_file, tmp_file))
        finally:
            os.remove(tmp_file)

    def test_multi_write_fasta(self):
        """
        End to end test for unique sequence in query
        """
        tmp_file = 'tmp_test_multi.output'
        expected_file = os.path.dirname(__file__) + '/multi.test.out'
        res = []
        for query in self.query:
            res.append(main.find_spacers(query, self.spacers))
        with open (tmp_file, 'w') as file_handle:
            main.write_fasta(res, file_handle)
        try:
            self.assertTrue(filecmp.cmp(expected_file, tmp_file))
        finally:
            os.remove(tmp_file)
