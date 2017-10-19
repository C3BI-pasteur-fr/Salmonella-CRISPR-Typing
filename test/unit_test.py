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

from salmonella_crispr_typing import main
from salmonella_crispr_typing.truncate_sequences import truncate_sequences
from salmonella_crispr_typing.extract_new import extract_new_spacers, _extract_sequences
from salmonella_crispr_typing.settings import LOCAL_DATA


class TestMain(unittest.TestCase):
    """
    Tests for main.py module
    """

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

    def test_list_spacers(self):
        """
        Test for listing of spacers
        """
        found_spacers = main.list_spacers(self.query, self.spacers)
        self.assertEqual(found_spacers[:6], "query1")

    def test_find_spacers(self):
        """
        Test for finding spacers in query with replacement
        """
        res = main.find_spacers(self.query[0], self.spacers)
        self.assertEqual(str(res.seq), "aa-spacer1&tcccggg")
        self.assertEqual(res.id, "query1")
        res = main.find_spacers(self.query[1], self.spacers)
        self.assertEqual(str(res.seq), "aa-spacer2&tcccggg")
        self.assertEqual(res.id, "query2")

    def test_clean_sequences(self):
        """
        Test function to clean sequences
        """
        res = [main.find_spacers(self.query[0], self.spacers)]
        cleaned_seq = main.clean_sequences(res)
        self.assertEqual(str(cleaned_seq[0].seq), "aa-spacer1tcccggg")

    def test_uniq_write_fasta(self):
        """
        End to end test for unique sequence in query
        """
        tmp_file = 'tmp_test_uniq.output'
        expected_file = os.path.dirname(__file__) + '/uniq.test.out'
        res = [(main.find_spacers(self.query[0], self.spacers))]
        with open(tmp_file, 'w') as file_handle:
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
        with open(tmp_file, 'w') as file_handle:
            main.write_fasta(res, file_handle)
        try:
            self.assertTrue(filecmp.cmp(expected_file, tmp_file))
        finally:
            os.remove(tmp_file)


class TestTruncateSequences(unittest.TestCase):
    """
    Tests for truncate_sequences.py module
    """

    def test_truncation(self):
        """
        End to end test for truncation of long genomic regions
        """
        tmp_file = 'tmp_test_trunc.output'
        input_file = os.path.dirname(__file__) + '/long_sequence.fa'
        expected_file = os.path.dirname(__file__) + '/trunc_example.test.output'
        # Open files
        with open(input_file, 'r') as file_handle:
            query = list(SeqIO.parse(file_handle, "fasta"))
        with open(LOCAL_DATA + "/spacers_Salmonella.fa", "r") as file_handle:
            spacers = list(SeqIO.parse(file_handle, "fasta"))
        res_query = [main.find_spacers(query[0], spacers)]
        res_query = truncate_sequences(res_query)
        with open(tmp_file, 'w') as file_handle:
            main.write_fasta(res_query, file_handle)
        try:
            self.assertTrue(filecmp.cmp(expected_file, tmp_file))
        finally:
            os.remove(tmp_file)

class TestExtractNew(unittest.TestCase):
    """
    Tests for extract_new.py module
    """

    def setUp(self):
        """
        Set up manual sequences for extraction
        """
        self.seqs = []
        self.seqs.append(SeqRecord(Seq("aaat-DR&-SPAC&-DR&ccggg"), id="no_new", name="query1"))
        self.seqs.append(SeqRecord(Seq("aaat-DR&actgatag-DR&gtcccggg"), id="new", name="query2"))

    def test_extract_new_spacers(self):
        """
        Test main method of extract_new module
        """
        new_seqs = extract_new_spacers(self.seqs)
        self.assertEqual(str(new_seqs[0].seq), "actgatag")

    def test__extract_sequences(self):
        """
        Test for _extract_sequences method
        """
        new_seq1 = _extract_sequences(self.seqs[0])
        self.assertFalse(new_seq1)
        new_seq2 = _extract_sequences(self.seqs[1])
        self.assertEqual(str(new_seq2[0].seq), "actgatag")
