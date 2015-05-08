#!/usr/bin/env python
"""
test_MetaPathways_filter_input.py

Created by Niels Hanson on 2015-05-08.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import unittest
import imp
MetaPathways_filter_input = imp.load_source('MetaPathways_filter_input', '../python_scripts/MetaPathways_filter_input.py')

class filterInputTestSet(unittest.TestCase):
    """
    Unit tests for MetaPathways_filter_input.py
    """
    def setUp(self):
        pass
    def testIsAminoAcidSequence(self):
        """
        Unit test for the function testIsAminoAcidSequence
        """
        protein_seqs = ["ABCDEFGHIJKLMNOP",\
                        "ARNDCQEGHILKMFPS"]
        nucleotide_seqs = ["AAAAAAAAAAAAAAAA",\
                           "ATCGATCGCGCGCGCG",
                           "CATTGCCAATATTTAGTGCATCTTCAATATTTTTATATTATTTTGCAGTCTAAAGACAGCTAACTATACTAAATTATTCGACAG"]
        
        ## Test sequences
        for pro_seq in protein_seqs:
            self.assertTrue(MetaPathways_filter_input.isAminoAcidSequence(pro_seq))
        
        for nuc_seq in nucleotide_seqs:
            self.assertFalse(MetaPathways_filter_input.isAminoAcidSequence(nuc_seq))
    
    def testFilterSequence(self):
        """
        Unit test for filter_sequence
        """
        
        test_seqs = ["CATTGCCAATATTTAGTGCATCTTCAATATTTTTATATTATTTTGCAGTCTAAAGACAGCTAACTATACTAAATTATTCGACAG",\
                     "CATTGCCAATATTTAGTGCATCTTCAANNNNNNNNNGTCTAAAGACAGCTAACTATACTAAATTATTCGACAG",\
                     "CATTGCCAATATTTZGTGCATCTTCAATATTTTTATATTATTTTGCAGTCTAAAGACAGCTAACTATACTAAATTATTCGACAG"]
        exp_seqs = ["CATTGCCAATATTTAGTGCATCTTCAATATTTTTATATTATTTTGCAGTCTAAAGACAGCTAACTATACTAAATTATTCGACAG",\
                    "GTCTAAAGACAGCTAACTATACTAAATTATTCGACAG",\
                    "GTGCATCTTCAATATTTTTATATTATTTTGCAGTCTAAAGACAGCTAACTATACTAAATTATTCGACAG"]
        
        for i in range(len(test_seqs)):
            self.assertEqual(MetaPathways_filter_input.filter_sequence(test_seqs[i]), exp_seqs[i], \
            "Did not correctly filter")
    


if __name__ == '__main__':
    unittest.main()