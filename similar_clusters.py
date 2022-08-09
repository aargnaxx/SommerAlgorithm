from sys import maxsize
from LengthCounting import LengthCounting, dirs
from difference_function import difference_function

debug = True


def find_most_similar_within_more_frequent_clusters(sequence, clusters):
    difference, min_difference, index_of_most_similar = maxsize, maxsize, -1
    for i, cluster in enumerate(clusters):
        difference = difference_function(sequence, cluster)
        if difference < min_difference:
            min_difference = difference
            index_of_most_similar = i
    return index_of_most_similar, min_difference


def _test_find_most_similar_within_more_frequent_clusters():
    c = LengthCounting(dirs[0], True)
    file = c.select_biggest_file()
    clusters = list(c.count_reading_occurrences(file, 201).keys())
    result = find_most_similar_within_more_frequent_clusters(clusters[-1], clusters[:-1])
    print(f'index of most similar: {result[0]}, difference: {result[1]}, self index: {len(clusters)}')


if __name__ == '__main__':
    _test_find_most_similar_within_more_frequent_clusters()
