import numpy as np
from collections import Counter
from FileLoader import scan_dir, FileLoader
from chimera import is_variant_a_chimera
from similar_clusters import find_most_similar_within_more_frequent_clusters

DEBUG = False


class AmpliconClassification:
    def __init__(self, readings, source=None):
        self.step, self.highest_artifact_frequency = 0, 0
        self.alleles, self.artifacts, self.unclassified, self.to_process = {}, {}, {}, {}
        self.artifacts_big_diff, self.artifacts_small_diff = {}, {}
        self.chimeras, self.temp_small_diff, self.temp_big_diff, self.temp_other, self.temp_same = [], [], [], [], []
        self.num_readings = len(readings)
        self.frequencies = \
            dict(sorted(Counter([reading.seq for reading in readings]).items(), key=lambda x: -x[1]))
        self.source = source

    def step_one(self, singleton_frequency=(1 / 10)):
        for sequence in self.frequencies.keys():
            if self.frequencies[sequence] < (self.num_readings * singleton_frequency):
                self.artifacts.update({sequence: self.frequencies[sequence]})
            else:
                self.to_process.update({sequence: self.frequencies[sequence]})
        if DEBUG:
            # print(f'to process: {list(self.to_process)}')
            # print(f'artifacts: {list(self.artifacts)}')
            pass
        if len(self.to_process) < 1:
            raise IOError(f'Input file {self.source} has no readings of desired length')
        first = list(self.to_process.keys())[0]
        self.alleles.update({first: self.frequencies[first]})
        for i, sequence in enumerate(list(self.to_process)[2:]):
            if is_variant_a_chimera(sequence, list(self.to_process)[:i]):
                self.chimeras.append(sequence)
            else:
                index, diff = find_most_similar_within_more_frequent_clusters(sequence, list(self.frequencies)[:i])
                if diff < 3:
                    self.temp_small_diff.append(sequence)
                else:
                    self.temp_big_diff.append(sequence)
        if len(self.frequencies) < 2:
            print(f'len {len(self.frequencies)} source {self.source}')
            print(f'fre {self.frequencies}')
            print(f'all {self.alleles}')
            print(f'art {self.artifacts}')
        del self.to_process

    def step_two_small_diff(self, replicate):
        further_processing = []
        for sequence in self.temp_small_diff:
            if sequence not in replicate.frequencies.keys():
                self.artifacts_small_diff.update({sequence: self.frequencies[sequence]})
                self.artifacts.update({sequence: self.frequencies[sequence]})
            else:
                further_processing.append(sequence)
        self.temp_small_diff = further_processing

    def step_two_big_diff(self, replicate, others):
        same, other = [], []
        for sequence in self.temp_big_diff:
            if sequence in replicate.frequencies.keys():
                same.append(sequence)
            else:
                for individual in others:
                    if sequence in individual.frequencies.keys():
                        other.append(sequence)
                        break
                self.artifacts.update({sequence: self.frequencies[sequence]})
                self.artifacts_big_diff.update({sequence: self.frequencies[sequence]})
        del self.temp_big_diff
        self.temp_same = same
        self.temp_other = other

    def step_two_chimera(self, replicate):
        chimeras_for_further_processing = []
        for sequence in self.chimeras:
            if sequence in replicate.frequencies.keys() and sequence not in replicate.chimeras:
                chimeras_for_further_processing.append(sequence)
            else:
                self.artifacts.update({sequence: self.frequencies[sequence]})
        self.chimeras = chimeras_for_further_processing

    def calculate_highest_artifacts_frequency(self):
        self.artifacts = dict(sorted(self.artifacts.items(), key=lambda x: -x[1]))
        if len(self.artifacts) > 0:
            self.highest_artifact_frequency = list(self.artifacts.items())[0][1]
        else:
            self.highest_artifact_frequency = 0
        if DEBUG:
            print(f'Highest artifact frequency: {self.highest_artifact_frequency}/{self.num_readings}')

    def step_three_small_diff(self):
        for sequence in self.temp_small_diff:
            item = {sequence: self.frequencies[sequence]}
            if self.frequencies[sequence] < self.highest_artifact_frequency:
                self.unclassified.update(item)
            else:
                self.alleles.update(item)
        del self.temp_small_diff

    def step_three_big_diff(self, others):
        for sequence in self.temp_same:
            item = {sequence: self.frequencies[sequence]}
            if self.frequencies[sequence] < self.highest_artifact_frequency:
                self.unclassified.update(item)
            else:
                self.alleles.update(item)
        del self.temp_same
        for sequence in self.temp_other:
            item = {sequence: self.frequencies[sequence]}
            for individual in others:
                if sequence in individual.alleles.keys() or sequence in individual.unclassified.keys():
                    self.unclassified.update(item)
                else:
                    self.artifacts.update(item)
                    self.artifacts_big_diff.update({sequence: self.frequencies[sequence]})
        del self.temp_other

    def step_three_chimera(self, others):
        for sequence in self.chimeras:
            for individual in others:
                if sequence in individual.alleles.keys():
                    self.alleles.update({sequence: self.frequencies[sequence]})

    def step_three_no_chimera_no_replica(self, others):
        for sequence in self.temp_small_diff + self.temp_big_diff:
            allele_flag = False
            item = {sequence: self.frequencies[sequence]}
            for individual in others:
                if sequence in individual.alleles.keys():
                    allele_flag = True
                    break
            frequency_flag = self.frequencies[sequence] < self.highest_artifact_frequency
            if allele_flag and not frequency_flag:
                self.alleles.update(item)
            elif frequency_flag:
                self.artifacts.update(item)
            else:
                self.unclassified.update(item)

    def step_three_chimera_no_replica(self, others):
        for sequence in self.chimeras:
            found = False
            for individual in others:
                if sequence in individual.chimeras:
                    self.artifacts.update({sequence: self.frequencies[sequence]})
                    found = True
                    break
                elif sequence in individual.alleles.keys() or sequence in individual.unclassified.keys():
                    self.unclassified.update({sequence: self.frequencies[sequence]})
                    found = True
                    break
            if not found:
                self.artifacts.update({sequence: self.frequencies[sequence]})

    def get_alleles(self):
        return self.alleles

    def get_artifacts(self):
        return self.artifacts, self.artifacts_small_diff, self.artifacts_big_diff

    def get_unclassified(self):
        return self.unclassified


def loop_function_call(nested_list, function):
    import warnings
    warnings.filterwarnings('ignore')
    if DEBUG:
        print(f'looping {function}')
    for i in nested_list:
        nested_np = np.array(nested_list)
        # print(nested_np)
        others = [j[0] for j in np.delete(nested_np, np.where(nested_np == i, 1, 0))]
        # print(others)
        if len(i) > 1:
            if function in ['step_two_small_diff', 'step_two_chimera']:
                getattr(i[0], function)(i[1])
            elif function == 'step_two_big_diff':
                getattr(i[0], function)(i[1], others)
            elif function in ['step_three_big_diff', 'step_three_chimera']:
                getattr(i[0], function)(others)
            elif function == 'step_one':
                getattr(i[0], function)()
        if len(i) == 1:
            if function in ['step_three_no_chimera_no_replica', 'step_three_chimera_no_replica']:
                getattr(i[0], function)(others)
            elif function == 'step_one':
                getattr(i[0], function)()


def run_classification(directory, reading_length=201):
    files = scan_dir(directory)
    individuals = list()
    for file in files:
        o, r, f, p = get_original_and_replica(file.path, reading_length, True)
        o = AmpliconClassification(o, source=str(f))
        r = AmpliconClassification(r, source=str(p)) if r else None
        # if r:
        #     if o.num_readings < r.num_readings:
        #         temp = o
        #         o = r
        #         r = temp
        if o.num_readings > 50:
            individuals.append([o, r] if r else [o])
    loop_function_call(individuals, 'step_one')
    loop_function_call(individuals, 'step_two_small_diff')
    loop_function_call(individuals, 'step_two_big_diff')
    loop_function_call(individuals, 'step_two_chimera')
    loop_function_call(individuals, 'calculate_highest_artifacts_frequency')
    loop_function_call(individuals, 'step_three_small_diff')
    loop_function_call(individuals, 'step_three_big_diff')
    loop_function_call(individuals, 'step_three_chimera')
    loop_function_call(individuals, 'step_three_no_chimera_no_replica')
    loop_function_call(individuals, 'step_three_chimera_no_replica')
    return individuals


def get_original_and_replica(filename, length=201, get_replica_path=False):
    # !!! WARNING !!! path variable depending on folder structure
    replica_path = filename.replace('A', 'B', 2)
    if DEBUG:
        print(filename)
        print(replica_path)
    original = FileLoader(filename).exact_length(length)
    try:
        replica = FileLoader(replica_path).exact_length(length)
    except FileNotFoundError:
        if DEBUG:
            print(f'{replica_path} not found')
        replica = None
    if get_replica_path:
        return original, replica, filename, replica_path
    return original, replica


if __name__ == '__main__':
    alleles_final, artifacts_final, unclassified_final = {}, {}, {}
    art1, art2 = {}, {}
    results = run_classification('C:/Users/S/diplomski/readings/DRB_A', 194)
    # results = run_classification('C:/Users/S/diplomski/readings/DQB_A', 201)
    more_than_3 = []
    more_than_2 = []
    for amplicon_pair in results:
        amplicon = amplicon_pair[0]
        alleles = amplicon.get_alleles()
        if len(alleles) > 0:
            print(f'{amplicon.source[-50:]}, {len(alleles)} alleles')
            for allele in alleles:
                print(f'   {alleles[allele]},{allele}')
                if allele not in alleles_final.keys():
                    alleles_final.update({allele: [[amplicon.source[-50:], alleles[allele], amplicon.num_readings]]})
                else:
                    alleles_final.update({allele: alleles_final[allele] + [
                        [amplicon.source[-50:], alleles[allele], amplicon.num_readings]
                    ]})
        if len(alleles) > 2:
            more_than_3.append(amplicon.source)
        if len(alleles) > 1:
            more_than_2.append(amplicon.source)
        artifacts, artifacts_from_small_diff, artifacts_from_big_diff = amplicon.get_artifacts()
        artifacts_final.update(artifacts)
        amp_unclassified = amplicon.get_unclassified()
        for unclassified in amp_unclassified:
            if unclassified not in unclassified_final.keys():
                unclassified_final.update(
                    {unclassified: [[amplicon.source[-50:], amp_unclassified[unclassified], amplicon.num_readings]]}
                )
            else:
                unclassified_final.update({unclassified: unclassified_final[unclassified] + [
                    [amplicon.source[-50:], amp_unclassified[unclassified], amplicon.num_readings]]
                                           })
    [print() for i in range(5)]
    print(f'Total alleles found {len(alleles_final)}')
    [print(i) for i in alleles_final.keys()]
    print()
    print('List of alleles and where they were found: ')
    for allele in alleles_final.keys():
        print(allele)
        for member in alleles_final[allele]:
            print('    ', member)
        print()
    print(f'Total artifacts found {len(artifacts_final)}')
    # print(artifacts)
    print(f'Total unclassified found {len(unclassified_final)}')
    [print(i) for i in unclassified_final.keys()]
    print()
    print('List of unclassified and where they were found: ')
    for unc in unclassified_final.keys():
        print(unc)
        for member in unclassified_final[unc]:
            print('    ', member)

    for i in more_than_3:
        print(i)
    print()
    print()
    print()
    for i in more_than_2:
        print(i)
