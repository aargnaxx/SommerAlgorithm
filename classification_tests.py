from AmpliconClassification import AmpliconClassification, loop_function_call
from main import get_original_and_replica

DEBUG = False


def test_step_one():
    orig, _ = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    a.step_one()
    print(len(a.alleles), len(a.artifacts), len(a.chimeras), len(a.temp_small_diff), len(a.temp_big_diff))


def test_step_two_a():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    for i in [a, b]:
        i.step_one()
    a.step_two_small_diff(b)


def test_step_two_b():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, b, c, d]:
        i.step_one()
    a.step_two_big_diff(b, [c, d])


def test_step_two_c():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, b, c, d]:
        i.step_one()
    a.step_two_chimera(b)


def test_step_three_a():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    for i in [a, b]:
        i.step_one()
    a.step_two_small_diff(b)
    a.calculate_highest_artifacts_frequency()
    a.step_three_small_diff()


def test_step_three_b():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, b, c, d]:
        i.step_one()
    a.step_two_big_diff(b, [c, d])
    a.calculate_highest_artifacts_frequency()
    a.step_three_big_diff([c, d])


def test_step_three_c():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, b, c, d]:
        i.step_one()
    a.step_two_chimera(b)
    b.step_two_chimera(a)
    a.step_three_chimera([c, d])


def test_step_three_d():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, c, d]:
        i.step_one()
    a.step_three_no_chimera_no_replica([c, d])


def test_step_three_e():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, c, d]:
        i.step_one()
    a.step_three_chimera_no_replica([c, d])


def test_artifacts_frequency():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)
    for i in [a, b, c, d]:
        i.step_one()
    a.step_two_small_diff(b)
    a.step_two_big_diff(b, [c, d])
    a.step_two_chimera(b)
    a.calculate_highest_artifacts_frequency()
    print(a.highest_artifact_frequency)


def general_test():
    orig, rep = get_original_and_replica('../readings/DQB_A/208_S208_L001_CS1A-URS_DQB.extendedFrags.fastq')
    a = AmpliconClassification(orig)
    b = AmpliconClassification(rep)
    orig, rep = get_original_and_replica('../readings/DQB_A/228_S228_L001_CS1A-URS_DQB.extendedFrags.fastq')
    c = AmpliconClassification(orig)
    d = AmpliconClassification(rep)

    list_with_all = [[a, b], [c, d]]
    loop_function_call(list_with_all, 'step_one')
    loop_function_call(list_with_all, 'step_two_small_diff')
    loop_function_call(list_with_all, 'step_two_big_diff')
    loop_function_call(list_with_all, 'step_two_chimera')
    loop_function_call(list_with_all, 'calculate_highest_artifacts_frequency')
    loop_function_call(list_with_all, 'step_three_small_diff')
    loop_function_call(list_with_all, 'step_three_big_diff')
    loop_function_call(list_with_all, 'step_three_chimera')


if __name__ == '__main__':
    test_step_one()
    test_step_two_a()
    test_step_two_b()
    test_step_two_c()
    test_step_three_a()
    test_step_three_b()
    test_step_three_c()
    test_step_three_d()
    test_step_three_e()
    general_test()
