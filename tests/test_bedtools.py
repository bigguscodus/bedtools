from bedtools import Bedtools
import unittest
import hashlib
bedtools = Bedtools()
class TestBedtools(unittest.TestCase):
    def multi_test(self, filename):
        h = hashlib.sha256()
        b = bytearray(128 * 1024)
        mv = memoryview(b)
        with open(filename, 'rb', buffering=0) as f:
            for n in iter(lambda: f.readinto(mv), 0):
                h.update(mv[:n])
        return h.hexdigest()
    def test_sort(self,expected = "bedtools_sorted.bed",observed = "sorted_lamina.bed"):
        bedtools._read_bed(first_bed="lamina.bed")
        bedtools.sort()
        self.assertEqual(self.multi_test(filename=expected), self.multi_test(filename=observed))
if __name__ == '__main__':
    unittest.main()