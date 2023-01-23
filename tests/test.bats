#!/usr/bin/env bats

load "${BATS_LIBDIR}/bats-support/load.bash"
load "${BATS_LIBDIR}bats-assert/load.bash"

setup() {
    TEST_DATA_DIR="$(dirname ${BATS_TEST_FILENAME})/data"
    TEST_TEMP_DIR="$(mktemp -d)"
}

teardown() {
    rm -rf "$TEST_TEMP_DIR"
}

@test "jf test data present" {
    jfdir_contents=$(ls ${TEST_DATA_DIR}/jf/*.jf | wc -l)
    assert_equal $jfdir_contents 4
}

@test "fasta test data present" {
    fastadir_contents=$(ls ${TEST_DATA_DIR}/fasta/*.fasta | wc -l)
    assert_equal $fastadir_contents 4
}

@test "fastq test data present" {
    fastqdir_contents=$(ls ${TEST_DATA_DIR}/fastq/*.fastq | wc -l)
    assert_equal $fastqdir_contents 4
}

@test "d2ssect works with fasta and jf input" {
    d2ssect -l ${TEST_DATA_DIR}/jf/*.jf -f ${TEST_DATA_DIR}/fasta/*.fasta 2>/dev/null  > ${TEST_TEMP_DIR}/test_matrix.txt
    matrix_diff=$(diff ${TEST_TEMP_DIR}/test_matrix.txt ${TEST_DATA_DIR}/matrix.txt)
    assert_equal ${matrix_diff} ""
}


@test "d2ssect works with fastq and jf input" {
    d2ssect -l ${TEST_DATA_DIR}/jf/*.jf -f ${TEST_DATA_DIR}/fastq/*.fastq 2>/dev/null  > ${TEST_TEMP_DIR}/test_matrix.txt

    matrix_diff=$(diff ${TEST_TEMP_DIR}/test_matrix.txt ${TEST_DATA_DIR}/matrix.txt)
    assert_equal ${matrix_diff} ""
}

