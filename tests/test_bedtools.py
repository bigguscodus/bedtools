from bedtools import Bedtools
import unittest
import hashlib
import os
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
    def test_sort(self,expected = "files_for_sort/bedtools_sorted.bed",observed = "sorted_lamina.bed"):
        bedtools.sort(bed_path="files_for_sort/lamina.bed")
        self.assertEqual(self.multi_test(filename=expected), self.multi_test(filename=observed))
    def test_merge(self, expected = "files_for_merge/bedtools_merge.bed", observed = "merged_lamina_m.bed"):
        bedtools.merge(bed_path="files_for_merge/lamina_m.bed")
        self.assertEqual(self.multi_test(filename=expected), self.multi_test(filename=observed))
    def test_getfasta(self,expected = "files_for_getfasta/bedtools_getfasta.fa", observed = "fasta.fa"):
        bedtools.getfasta(path_to_fasta="files_for_getfasta/sequence.fasta",
                          path_to_bed="files_for_getfasta/sample.bed")
        self.assertEqual(self.multi_test(filename=expected), self.multi_test(filename=observed))
    def test_intesect(self, expected = "files_for_intersect/bedtools_intersected.bed", observed = "intersected_bed_1_bed_2.bed"):
        bedtools.intersect("files_for_intersect/bed_1.bed", "files_for_intersect/bed_2.bed", wo=None)
        self.assertEqual(self.multi_test(filename=expected), self.multi_test(filename=observed))

if __name__ == '__main__':
    unittest.main()
