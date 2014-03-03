import unittest
import subprocess
import os
import regionanalysis
import string


class TestRA(unittest.TestCase):

    def setUp(self):
        self.module_dir = os.path.dirname(
            os.path.realpath(regionanalysis.__file__))
        self.example_dir = os.path.join(self.module_dir, "example/")
        self.exception_dir = os.path.join(self.module_dir, "exceptions/")
        print(self.module_dir)

    def test_with_header(self):
        input_example = os.path.join(self.example_dir, "test_with_header.bed")
        cmds = ["region_analysis.py", "-d", "ensembl",
                "-g", "mm10", "-i", input_example, "-r"]
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE, bufsize=1)
        stdout, stderr = p.communicate()
        test_output = open(input_example + ".annotated").read()
        test_ref = open(input_example + ".annotated.out").read()
        self.assertMultiLineEqual(test_output, test_ref)
        test_output = open(input_example + ".full.annotated").read()
        test_ref = open(input_example + ".full.annotated.out").read()
        self.assertMultiLineEqual(test_output, test_ref)
        test_output = open(input_example + ".full.annotated.json").read()
        test_ref = open(input_example + ".full.annotated.json.out").read()
        self.assertMultiLineEqual(test_output, test_ref)

    def test_without_header(self):
        input_example = os.path.join(
            self.example_dir, "test_without_header.bed")
        cmds = ["region_analysis.py", "-d", "ensembl",
                "-g", "mm10", "-i", input_example]
        p = subprocess.Popen(cmds, stdout=subprocess.PIPE, bufsize=1)
        stdout, stderr = p.communicate()
        test_output = open(input_example + ".annotated").read()
        test_ref = open(input_example + ".annotated.out").read()
        self.assertMultiLineEqual(test_output, test_ref)
        test_output = open(input_example + ".full.annotated").read()
        test_ref = open(input_example + ".full.annotated.out").read()
        self.assertMultiLineEqual(test_output, test_ref)
        test_output = open(input_example + ".full.annotated.json").read()
        test_ref = open(input_example + ".full.annotated.json.out").read()
        self.assertMultiLineEqual(test_output, test_ref)

    def test_exceptions(self):
        # in subprocess.Popen, shell=True, executable="bash" were used, because in UBUNTU, /bin/sh is not bash!
        # null input
        null_msg = "Please assign proper input file!"
        cmds = ["region_analysis.py"]
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            bufsize=1, universal_newlines=True, shell=False)
        stdout, stderr = p.communicate()
        # split stderr output into list
        stderr_list = string.split(stderr, sep="\n")
        # check if the expected error message is in the stderr output
        is_exception = any([(null_msg in y) for y in stderr_list])
        self.assertEqual(is_exception, True)
        del p, stdout, stderr
        # try different input files
        genomes = ["dumb", "mm10", "mm10"]
        msgs = ["dumb not in the genome database!",
                "Error in input file! Please check the format!",
                "Error in input file! Please check the format!"]
        for i in range(0, 3):
            input_example = os.path.join(self.exception_dir, "%d.bed" % i)
            cmds = ["region_analysis.py", "-d", "ensembl",
                    "-g", genomes[i], "-i", input_example]
            print(cmds)
            p = subprocess.Popen(
                cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                bufsize=1, universal_newlines=True, shell=False)
            stdout, stderr = p.communicate()
            self.assertEqual(msgs[i], string.strip(stderr))
            del p, stdout, stderr


if __name__ == '__main__':
    unittest.main()
