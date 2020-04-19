import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("-fb", "--first_bed", help="Path to the first bed file")
parser.add_argument("-sb", "--second_bed", help="path to the second bed file")
args = parser.parse_args()


def _merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    interval_index = 0
    # print(sorted_intervals)
    for i in sorted_intervals:
        if i[0] > sorted_intervals[interval_index][1]:
            interval_index += 1
            sorted_intervals[interval_index] = i
        else:
            sorted_intervals[interval_index] = [sorted_intervals[interval_index][0], i[1]]
    # print(sorted_intervals)
    return sorted_intervals[:interval_index + 1]

class Bedtools:
    def __init__(self):
        pass

    def _create_bed_generator(self, bed_path):
        if bed_path is not None:
            with open(bed_path) as file:
                self.bed_generator = (str(line).strip().split("\t") for line in file.readlines())
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

    def merge(self,bed_path=args.first_bed):
        reads_for_names = self._create_bed_generator(bed_path)
        reads_for_values = self._create_bed_generator(bed_path)
        chr_name_dict = OrderedDict((read[0],list()) for read in reads_for_names)
        for read in reads_for_values:
            chr_name_dict[read[0]].append([int(read[1]), int(read[2])])
        name = self._parse_bed_name(bed_path=bed_path)
        merged_name = f"merged_{name}"
        with open(merged_name,"w+") as out:
            for i, (key, value) in enumerate(chr_name_dict.items()):
                for z in _merge_intervals(intervals=value):
                    p = "\t".join([key, str(z[0]), str(z[1])])
                    out.write(f"{p}\n")
    def intersect(self):
        pass

    def getfasta(self):
        pass
)
