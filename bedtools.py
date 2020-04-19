import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-fb", "--first_bed", help="Path to the first bed file")
parser.add_argument("-sb", "--second_bed", help="path to the second bed file")
args = parser.parse_args()
class Bedtools:
    def __init__(self):
        pass
    def _read_bed(self, first_bed = args.first_bed, second_bed = args.second_bed):
        if (first_bed is not None):
            with open(first_bed) as file_1:
                self.first_bed_generator = (str(line).strip().split("\t") for line in file_1.readlines()[1:])
        else:
            self.first_bed_generator = None
        if (second_bed is not None):
            with open(second_bed) as file_2:
                self.second_bed_generator = (str(line).strip().split("\t") for line in file_2.readlines()[1:])
        else:
            self.second_bed_generator = None
        print(self.first_bed_generator)
        print(self.second_bed_generator)
        #print(self.second_bed_generator)
    def sort(bed_to_sort = args.first_bed):
        pass
    def substract(self):
        pass
    def merge(self):
        pass
    def intersect(self):
        pass
    def getfasta(self):
        pass
if __name__ == "__main__":
    bedtools = Bedtools()
    bedtools._read_bed()
