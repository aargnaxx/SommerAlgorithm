import os
from itertools import chain
from Bio import SeqIO

DEBUG = False


class FileLoader:
    def __init__(self, files, format_="fastq"):
        self.readings = list()
        if isinstance(files, list):
            self.files = files
            if DEBUG:
                print(self.files)
        else:
            self.files = [files]
        self.format = format_

    def filter_readings(self, f):
        self.readings = list(chain.from_iterable(
            [[rec for rec in SeqIO.parse(file, self.format)
              if f(rec)] for file in self.files]
        ))
        if DEBUG:
            print(f'read {len(self.files)} files')
            print(f'read {len(self.readings)} sequences')
        return self.readings

    def all(self):
        if DEBUG:
            print('all readings')
        return self.filter_readings((lambda x: True))

    def longer_than(self, filtered_length):
        if DEBUG:
            print(f'length longer than {filtered_length}')
        return self.filter_readings(lambda x: len(x.seq) > filtered_length)

    def exact_length(self, filtered_length):
        if DEBUG:
            print(f'exact length {filtered_length}')
        return self.filter_readings(lambda x: len(x.seq) == filtered_length)

    def length_between(self, filtered_length):
        if DEBUG:
            print(f'min {min(filtered_length) + 1}, max {max(filtered_length) + 1}')
        return self.filter_readings(lambda x: (min(filtered_length) < len(x.seq) < max(filtered_length)))


def scan_dir(path):
    files = [f for f in os.scandir(path)]
    return files


def get_file_size(file):
    return os.path.getsize(file)


def as_dict(sequences):
    return SeqIO.to_dict(sequences)


if __name__ == '__main__':
    file_ = 'C:/Users/S/diplomski/ocitanja/DQB_A/193_S193_L001_CS1A-URS_DQB.extendedFrags.fastq'
    filtered = FileLoader(file_)
    readings = filtered.all()
    print('example seq:')
    print(readings[0].seq)


# TODO filter file by reading length and quantity