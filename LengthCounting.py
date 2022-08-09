from collections import Counter
from FileLoader import FileLoader, scan_dir, get_file_size

dirs = [
    'C:/Users/S/diplomski/readings/DQB_A',
    'C:/Users/S/diplomski/readings/DQB_B',
    'C:/Users/S/diplomski/readings/DRB_A',
    'C:/Users/S/diplomski/readings/DRB_B',
]


class LengthCounting:
    def __init__(self, directory, debug=False):
        self.directory = directory
        self.debug = debug
        self.count = {}

    def count_sequence_lengths(self):
        files = scan_dir(self.directory)
        reader = FileLoader(files)
        count = Counter(
            [len(reading) for reading in reader.all()]
        )
        if self.debug:
            print(count)
        return count

    def count_reading_occurrences(self, file_, length_):
        reader = FileLoader(file_)
        readings = reader.exact_length(length_)
        count = Counter([reading.seq for reading in readings])
        self.count = count
        return self.count

    def display_count(self, display_first_=20):
        for reading in sorted(list(self.count.keys()), key=lambda x: -self.count[x])[:display_first_]:
            print(f'{self.count[reading]:{5}},\t {reading}')

    def select_biggest_file(self):
        files = scan_dir(self.directory)
        max_size, index = 0, 0
        filename = ''
        for i, f in enumerate(files):
            current_size = get_file_size(f)
            if current_size > max_size:
                max_size = current_size
                index = i
                filename = f
        if self.debug:
            print(f'largest file {filename} at index {index}')
        return files[index]


def _test_length_counting_():
    c = LengthCounting(dirs[0], True)
    file = c.select_biggest_file()
    c.count_reading_occurrences(file, 201)
    c.display_count()


def _test_sequence_lengths_():
    c = LengthCounting(dirs[2], True)
    c.count_sequence_lengths()


if __name__ == '__main__':
    # _test_length_counting_()
    _test_sequence_lengths_()
