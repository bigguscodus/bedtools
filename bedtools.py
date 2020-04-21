import argparse
import itertools
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("-fb", "--first_bed", help="Path to the first bed file")
parser.add_argument("-sb", "--second_bed", help="path to the second bed file")
parser.add_argument("-fi", "--fasta_file", help="path to fasta file")
parser.add_argument("-bed", "--bed_file", help="path to sorted bed file")
parser.add_argument("-wo", "--wo", help="Equal to bedtools intersect wo")
args = parser.parse_args()


class Bedtools:
    def __init__(self):
        pass

    def _intervals_extract(self,iterable):
        iterable = sorted(set(iterable))
        for key, group in itertools.groupby(enumerate(iterable),
                                            lambda t: t[1] - t[0]):
            group = list(group)
            yield [group[0][1], group[-1][1]]

    def _create_OrderDict(self,bed_generator_1, bed_generator_2):
        result = OrderedDict((read[0], list()) for read in bed_generator_1)
        for read in bed_generator_2:
            result[read[0]].extend([i for i in range(int(read[1]), int(read[2]) + 1)])
        return result

    def _intersect(self,first_bed_dict, second_bed_dict, wo):
        for ikey, ivalue in first_bed_dict.items():
            for jkey, jvalue in second_bed_dict.items():
                if ikey == jkey:
                    intersect = self._intervals_extract(value for value in ivalue if value in jvalue)
                    for interval in intersect:
                        if wo is None:
                            yield f"{ikey}\t{interval[0]}\t{interval[1]}"
                        else:
                            yield f"{ikey}\t{ivalue[0]}\t{ivalue[-1]}\t{jkey}\t{jvalue[0]}\t{jvalue[1]}\t{interval[1] - interval[0]}"
    def _merge_intervals(self, intervals):
        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        interval_index = 0
        for i in sorted_intervals:
            if i[0] > sorted_intervals[interval_index][1]:
                interval_index += 1
                sorted_intervals[interval_index] = i
            else:
                sorted_intervals[interval_index] = [sorted_intervals[interval_index][0], i[1]]
        return sorted_intervals[:interval_index + 1]

    def _read_fasta(self, fasta):
        name, seq = None, []
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

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

    def substract(self, first_bed_path=args.first_bed, second_bed_path=args.second_bed):
        bed_generator_1 = self._create_bed_generator(first_bed_path)
        bed_generator_1_2 = self._create_bed_generator(first_bed_path)
        bed_generator_2 = self._create_bed_generator(second_bed_path)
        bed_generator_2_2 = self._create_bed_generator(second_bed_path)
        first_bed_dict = OrderedDict((read[0], list()) for read in bed_generator_1)
        for read in bed_generator_1_2:
            first_bed_dict[read[0]].extend([i for i in range(int(read[1]), int(read[2]) + 1)])
        second_bed_dict = OrderedDict((read[0], list()) for read in bed_generator_2)
        for read in bed_generator_2_2:
            second_bed_dict[read[0]].extend([i for i in range(int(read[1]), int(read[2]) + 1)])
        used_keys = list()
        with open("substract.bed", "w+") as out:
            for i, (key, value) in enumerate(first_bed_dict.items()):
                for j, (jkey, jvalue) in enumerate(second_bed_dict.items()):
                    if key == jkey:
                        if value[0] < jvalue[0]:
                            substract = [a + 1 for a in value if a not in jvalue]
                            out.write(f"{key}\t{substract[0] - 1}\t{substract[-1]}\n")
                        if value[0] >= jvalue[0]:
                            substract = [a for a in value if a not in jvalue]
                            out.write(f"{key}\t{substract[0] - 1}\t{substract[-1]}\n")
                    if key not in second_bed_dict and key not in used_keys:
                        out.write(f"{key}\t{value[0]}\t{value[-1]}\n")
                        used_keys.append(key)

    def merge(self, bed_path=args.first_bed):
        reads_for_names, reads_for_values = self._create_bed_generator(bed_path), self._create_bed_generator(bed_path)
        chr_bed_dict = OrderedDict((read[0], list()) for read in reads_for_names)
        for read in reads_for_values:
            chr_bed_dict[read[0]].append([int(read[1]), int(read[2])])
        name = self._parse_bed_name(bed_path=bed_path)
        merged_name = f"merged_{name}"
        with open(merged_name, "w+") as out:
            for i, (key, value) in enumerate(chr_bed_dict.items()):
                for interval in self._merge_intervals(intervals=value):
                    line = "\t".join([key, str(interval[0]), str(interval[1])])
                    out.write(f"{line}\n")

    def intersect(self, first_bed_path=args.first_bed, second_bed_path=args.second_bed, wo=args.wo):
        bed_generator_1 = self._create_bed_generator(first_bed_path)
        bed_generator_2 = self._create_bed_generator(second_bed_path)
        bed_generator_1_2 = self._create_bed_generator(first_bed_path)
        bed_generator_2_2 = self._create_bed_generator(second_bed_path)
        first_bed_dict = self._create_OrderDict(bed_generator_1, bed_generator_1_2)
        second_bed_dict = self._create_OrderDict(bed_generator_2, bed_generator_2_2)
        part_name_1 = self._parse_bed_name(first_bed_path)
        part_name_2 = self._parse_bed_name(second_bed_path)
        with open(f"intersected_{part_name_1.rstrip(part_name_1[-4:])}_{part_name_2.rstrip(part_name_2[-4:])}.bed","w+") as file:
            lines = self._intersect(first_bed_dict=first_bed_dict,second_bed_dict=second_bed_dict,wo=wo)
            for line in lines:
                file.write(f"{line}\n")

    def getfasta(self, path_to_fasta=args.fasta_file, path_to_bed=args.bed_file):
        bed_generator = self._create_bed_generator(bed_path=path_to_bed)
        with open(path_to_fasta) as fasta:
            for name, seq in self._read_fasta(fasta):
                for read in bed_generator:
                    if read[0] in name:
                        line = [f">{read[0]}:{read[1]}-{read[2]}", seq[int(read[1]):int(read[2])]]
                        str_line = "\n".join(line)
                        with open('fasta.fa', "a+") as out:
                            out.write(f"{str_line}\n")

