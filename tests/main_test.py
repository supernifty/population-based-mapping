
import collections
import StringIO
import sys
import unittest

sys.path.append('./src/')
import main

class TestMain(unittest.TestCase):

    def test_write_tsv(self):
        u = collections.defaultdict(int)
        m = collections.defaultdict(int)
        o = collections.defaultdict(int)
        u[1] = 0.1
        u[10] = 0.2
        m[1] = 0.3
        m[10] = 0.4
        o[3] = 0.5
        out = StringIO.StringIO()
        main.write_tsv(out, u, m, o)
        lines = out.getvalue().split('\n')
        # ['Position\tUnique\tMulti\tOther', '1\t0.100\t0.300\t0.000', '3\t0.000\t0.000\t0.500', '10\t0.200\t0.400\t0.000', '']
        self.assertEqual(5, len(lines))
        self.assertEqual('1\t0.100\t0.300\t0.000', lines[1])
        self.assertEqual('3\t0.000\t0.000\t0.500', lines[2])
        self.assertEqual('10\t0.200\t0.400\t0.000', lines[3])
        self.assertEqual('', lines[4])

    def test_write_stats(self):
        u = collections.defaultdict(int)
        m = collections.defaultdict(int)
        o = collections.defaultdict(int)
        u[1] = 0.1
        u[10] = 0.2
        m[1] = 0.3
        m[10] = 0.4
        o[3] = 0.5
        out = StringIO.StringIO()
        main.write_stats(out, u, m, o)
        lines = out.getvalue().split('\n')
        # ['Type\tCount\tMax\tMin\tMean', 'Unique\t2\t0.2\t0.1\t0.15', 'Multi\t2\t0.4\t0.3\t0.35', 'Other\t1\t0.5\t0.5\t0.5', '']
        self.assertEqual(5, len(lines))
        self.assertEqual('Unique\t2\t0.2\t0.1\t0.15', lines[1])
        self.assertEqual('Multi\t2\t0.4\t0.3\t0.35', lines[2])
        self.assertEqual('Other\t1\t0.5\t0.5\t0.5', lines[3])
        self.assertEqual('', lines[4])

if __name__ == '__main__':
    unittest.main()
