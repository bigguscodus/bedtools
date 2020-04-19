import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-fb", "--first_bed", help="Path to the first bed file")
parser.add_argument("-sb", "--second_bed", help="path to the second bed file")
args = parser.parse_args()


class Bedtools:
    def __init__(self):
        pass

    def _create_bed_generator(self, bed_path):
        if bed_path is not None:
            with open(bed_path) as file:
                self.bed_generator = (str(line).strip().split("\t") for line in file.readlines()[1:])
        else:
            self.bed_generator = None
        return self.bed_generator

    def _parse_bed_name(self, bed_path: str):
        if "/" in bed_path:
            name = bed_path.split('/')[-1]
        else:
            name = bed_path
        return name

    def sort(self, bed_path=args.first_bed):
        bed_generator = self._create_bed_generator(bed_path=bed_path)
        if bed_generator is None:
            raise ValueError("Need .bed file for sorting")
        name = self._parse_bed_name(bed_path=bed_path)
        sorted_name = f"sorted_{name}"
        with open(sorted_name, "w+") as sorted_out:
            for read in sorted(bed_generator, key=lambda x: (x[0], int(x[1]))):
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

