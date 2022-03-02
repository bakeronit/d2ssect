import unittest
import tempfile
import filecmp
import os
import d2ssect.jellyfish as jellyfish

## import other functions to be tested

class TestJellyfish(unittest.TestCase):

        
    def test_jfcount(self):
        inputf="data/fasta/DI-1-1_S6.fasta"
        outputf=tempfile.mkstemp()[1]

        print(outputf)
        os.system("cat "+outputf)
        jellyfish.count(inputf,outputf)

        expected_output="data/jf/DI-1-1_S6.jf"

        # Note files are compared on size 
        # only as jellyfish will not produce identical files from separate runs
        self.assertTrue(os.path.getsize(outputf),os.path.getsize(expected_output))

        
if __name__ == '__main__':
    unittest.main()
