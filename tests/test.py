import unittest
import tempfile
import filecmp
import os, shutil
import d2ssect.jellyfish as jellyfish

## import other functions to be tested

class TestJellyfish(unittest.TestCase):

    def setUp(self) -> None:
        self.test_dir = tempfile.mkdtemp()
        return super().setUp()
        
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

        
    def test_jfdump(self):
        inputf="data/jf/D1-1-1_S6.jf"
        outputf=os.path.join(self.test_dir,"D1-1-1_S6.nkz")
        
        jellyfish.dump(inputf, outputf)
        expected_output="data/nkz/DI-1-1_S6.nkz"

        self.assertTrue(filecmp.cmp(outputf, expected_output))


    def tearDown(self) -> None:
        shutil.rmtree(self.test_dir)
        return super().tearDown()

if __name__ == '__main__':
    unittest.main()
