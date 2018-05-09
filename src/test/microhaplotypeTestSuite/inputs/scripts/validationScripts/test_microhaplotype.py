import pytest
import pysam
import shlex
import subprocess
import filecmp

@pytest.fixture
def inputs_path():
    return './test_microhaplotype/input_files'


@pytest.fixture
def output_path():
    return './test_microhaplotype/output_files'


@pytest.fixture
def default_obj(inputs_path, output_path):
    from main import Microhaplotyper
    return Microhaplotyper(bam=inputs_path + '/test_coverage_10_all_reads_cover.bam',
                           bed=inputs_path + '/test_default.bed',
                           info=inputs_path + '/test_info.txt',
                           out=output_path + '/test_output.txt',
                           analytical_threshold=0.01)
@pytest.fixture
def default_obj_no_threshold(inputs_path, output_path):
    from main import Microhaplotyper
    return Microhaplotyper(bam=inputs_path + '/test_coverage_10_all_reads_cover.bam',
                           bed=inputs_path + '/test_default.bed',
                           info=inputs_path + '/test_info.txt',
                           out=output_path + '/test_output.txt')

@pytest.fixture
def bed_obj(inputs_path):
    from bed import Bed
    return Bed(bed_file_path=inputs_path + '/test_default.bed')

@pytest.fixture
def bad_bed_obj(inputs_path):
    from bed import Bed
    return Bed(bed_file_path=inputs_path + '/test_bad_bed_less_cols.bed')

@pytest.fixture
def info_obj(inputs_path):
    from info import Info
    return Info(info_file_path=inputs_path + '/test_info.txt')

@pytest.fixture
def microhap_obj():
    from microhaplotype import Microhaplotype
    return Microhaplotype(name='mitotrial-001')

@pytest.fixture
def microhap_processed_obj(microhap_obj, inputs_path):
    bam_file_obj = pysam.AlignmentFile(inputs_path + '/test_coverage_10_all_reads_cover.bam', 'rb')
    microhap_obj.process_overlapping_reads('chrM', 149, bam_file_obj)
    microhap_obj.process_overlapping_reads('chrM', 199, bam_file_obj)
    microhap_obj.process_overlapping_reads('chrM', 249, bam_file_obj)
    return microhap_obj

@pytest.fixture
def filter_obj_default():
    from filters import Filters
    return Filters(analytical_threshold=0.02, abs_analytical_threshold=10, min_mapq=60, max_mapq=254,
                   remove_deletions=True)

@pytest.fixture
def filter_obj_alt():
    from filters import Filters
    return Filters(analytical_threshold=0.02, abs_analytical_threshold=11, min_mapq=60, max_mapq=254,
                   remove_deletions=True)

@pytest.fixture
def argument_parsing_base(inputs_path, output_path):
    return 'python main.py' + \
                   ' --bam ' + inputs_path + '/test_coverage_10_all_reads_cover.bam' + \
                   ' --bed ' + inputs_path + '/test_default.bed'

@pytest.fixture
def allele_representation_test():
    from allele import Allele
    return Allele(sequence='ATGC', coverage_total=10, coverage_forward=5, coverage_reverse=5, total_mapping_quality=700)

@pytest.fixture
def output_obj(default_obj_no_threshold):
    default_obj_no_threshold.get_microhaplotypes()
    from output_microhaplotypes import OutputMicrohaplotypes
    return OutputMicrohaplotypes(default_obj_no_threshold)

# ========================Argument Parser=====================
def test_argument_sanity(inputs_path, argument_parsing_base):
    command_to_execute = str(argument_parsing_base)
    assert subprocess.check_call(shlex.split(command_to_execute)) == 0

def test_argument_sanity2(argument_parsing_base, inputs_path, output_path):
    command_to_execute = str(argument_parsing_base) + \
                         ' --info ' + str(inputs_path) + '/test_info.txt' + \
                         ' --out ' + str(output_path) + '/test_output.txt' + \
                         ' --mincov 0.03'
    assert subprocess.check_call(shlex.split(command_to_execute)) == 0


def test_argument_mincov_allowed(argument_parsing_base):
    for mincov_value in [0.0, 1.0, 0.5, 0.00006]:
        command_to_execute = str(argument_parsing_base) + \
                             ' --mincov ' + str(mincov_value)
        assert subprocess.check_call(shlex.split(command_to_execute)) == 0

def test_argument_mincov_not_allowed(argument_parsing_base):
    for mincov_value in [-1, 10, 3]:
        command_to_execute = str(argument_parsing_base) + \
                             ' --mincov ' + str(mincov_value)
        with pytest.raises(subprocess.CalledProcessError):
            subprocess.check_call(shlex.split(command_to_execute))


# ========================Filter Tests========================
def test_filters(filter_obj_default):
    assert filter_obj_default.check_analytical_threshold(100, 1) is False
    assert filter_obj_default.check_analytical_threshold(100, 2) is True
    assert filter_obj_default.check_absolute_analytical_threshold(9) is False
    assert filter_obj_default.check_absolute_analytical_threshold(10) is True
    assert filter_obj_default.check_min_mapq(59) is False
    assert filter_obj_default.check_min_mapq(60) is True
    assert filter_obj_default.check_max_mapq(255) is False
    assert filter_obj_default.check_max_mapq(254) is True
    assert filter_obj_default.check_deletions('A_G') is False
    assert filter_obj_default.check_deletions('ATG') is True
    assert filter_obj_default.get_analytical_threshold(100) == 2

# ========================Microhaplotype Tests========================

def test_process_overlapping_reads(microhap_obj, inputs_path):
    bam_file_obj = pysam.AlignmentFile(inputs_path + '/test_coverage_10_all_reads_cover.bam', 'rb')
    microhap_obj.process_overlapping_reads('chrM', 149, bam_file_obj)
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].mh_sequence == 'G'
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].is_reverse is False
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].mapping_quality == 73
    microhap_obj.process_overlapping_reads('chrM', 199, bam_file_obj)
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].mh_sequence == 'GT'
    microhap_obj.process_overlapping_reads('chrM', 249, bam_file_obj)
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].mh_sequence == 'GTG'
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].is_reverse is False
    assert microhap_obj.reads['contributor1_chrM_101_1_0_0_0_0_0:3:0_0:0:0_2'].mapping_quality == 73


def test_get_alleles(microhap_processed_obj):
    total_reads, short_reads = microhap_processed_obj.get_alleles(3)
    assert microhap_processed_obj.alleles['GTG'].sequence == 'GTG'
    assert microhap_processed_obj.alleles['GTG'].coverage_total == 10
    assert microhap_processed_obj.alleles['GTG'].coverage_forward == 7
    assert microhap_processed_obj.alleles['GTG'].coverage_reverse == 3
    assert microhap_processed_obj.alleles['GTG'].total_mapping_quality == 701
    assert total_reads == 10
    assert short_reads == 0


def test_calculate_avg_mapping_quality(microhap_processed_obj):
    microhap_processed_obj.get_alleles(3)
    microhap_processed_obj.calculate_avg_mapping_quality()
    assert microhap_processed_obj.alleles['GTG'].avg_mapping_quality == 70


def test_calculate_total_coverage(microhap_processed_obj):
    microhap_processed_obj.get_alleles(3)
    microhap_processed_obj.calculate_total_coverage()
    assert microhap_processed_obj.total_coverage == 10

def test_filter_alleles(microhap_processed_obj, filter_obj_default, filter_obj_alt):
    microhap_processed_obj.get_alleles(3)
    microhap_processed_obj.calculate_avg_mapping_quality()
    microhap_processed_obj.calculate_total_coverage()
    microhap_processed_obj.filter_alleles_and_set_analytical_threshold(filter_obj_alt)
    with pytest.raises(KeyError):
        test = microhap_processed_obj.alleles_filtered['GTG']
    microhap_processed_obj.filter_alleles_and_set_analytical_threshold(filter_obj_default)
    assert microhap_processed_obj.alleles_filtered['GTG'].sequence == 'GTG'
    assert microhap_processed_obj.alleles_filtered['GTG'].coverage_total == 10
    assert microhap_processed_obj.alleles_filtered['GTG'].coverage_forward == 7
    assert microhap_processed_obj.alleles_filtered['GTG'].coverage_reverse == 3
    assert microhap_processed_obj.alleles_filtered['GTG'].total_mapping_quality == 701
    assert microhap_processed_obj.alleles_filtered['GTG'].avg_mapping_quality == 70
    assert microhap_processed_obj.analytical_threshold == 0

def test_get_allele_coverage_str(microhap_processed_obj, filter_obj_default):
    microhap_processed_obj.get_alleles(3)
    microhap_processed_obj.calculate_avg_mapping_quality()
    microhap_processed_obj.calculate_total_coverage()
    microhap_processed_obj.filter_alleles_and_set_analytical_threshold(filter_obj_default)
    assert str(microhap_processed_obj.get_allele_coverage_str()) == '[GTG]10'
    assert str(microhap_processed_obj.get_allele_coverage_str_detailed()) == '[GTG]10(3:7:70)'



# ========================Microhaplotyper Tests=======================

def test_default_analytical_threshold(default_obj_no_threshold):
    """ Tests if the default value for analytical threshold """
    assert default_obj_no_threshold.analytical_threshold == 0.02


def test_mh_output_concordance(output_obj, inputs_path, output_path):
    output_obj.output_microhaplotypes()
    assert filecmp.cmp(inputs_path + '/test_output.txt', output_path + '/test_output.txt')


# =========================Bed Tests ================================
def test_bed_header_filter(bed_obj):
    """ Tests ignoring of header lines """
    assert bed_obj._bed_header_line_check("track something") == False
    assert bed_obj._bed_header_line_check("# something") == False

def test_get_bed_dict(bed_obj):
    bed_obj.get_bed_dict()
    assert bed_obj.bed_dict['mitotrial-001'] == [('chrM', 149), ('chrM', 199), ('chrM', 249)]

def test_bad_bed(bad_bed_obj):
    with pytest.raises(SystemExit):
        bad_bed_obj.get_bed_dict()

# ========================Info File Tests ===========================
def test_get_info_dict(info_obj):
    """ Test info file dict """
    info_obj.get_info_dict()
    assert info_obj.info_dict['mitotrial-001']['Ae'] == 5.54822
    assert info_obj.info_dict['mitotrial-001']['ProbOfDetMix'] == 0.37857
    assert info_obj.info_dict['mitotrial-001']['RSIDList'] == 'rs1927847/rs9536429'

def test_add_info(info_obj, microhap_obj):
    info_obj.get_info_dict()
    info_obj.add_info(microhap_obj)
    assert microhap_obj.ae == 5.54822
    assert microhap_obj.prob_det_mix == 0.37857
    assert microhap_obj.snps_rsid_list == 'rs1927847/rs9536429'

# =======================Allele Tests===============================
def test_allele_class(allele_representation_test):
    allele_representation_test.compute_avg_mapping_quality()
    assert allele_representation_test.avg_mapping_quality == 70
    assert str(allele_representation_test.get_detailed_representation()) == '[ATGC]10(5:5:70)'
    assert str(allele_representation_test.get_simple_representation()) == '[ATGC]10'
    allele_representation_test.coverage_total = 0
    with pytest.raises(SystemExit):
        allele_representation_test.compute_avg_mapping_quality()
