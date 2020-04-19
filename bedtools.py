import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-fb", "--first_bed", help="Path to the first bed file")
parser.add_argument("-sb", "--second_bed", help="path to the second bed file")
args = parser.parse_args()
class Bedtools:
    def __init__(self):
        pass
    def _read_bed(self, first_bed = args.first_bed, second_bed = args.second_bed):
        self.first_path = first_bed
        self.second_path = second_bed
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
    def sort(self):
        if (self.first_bed_generator is None):
            raise ValueError("Need .bed file for sorting")
        name = f"sorted_{self.first_path.split('/')[-1]}"
        with open(name,"w+") as sorted_out:
            for read in sorted(self.first_bed_generator, key=lambda x: (x[0], int(x[1]))):
                read = "\t".join(read)
                sorted_out.write(f"{read}\n")
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
    bedtools.sort()
