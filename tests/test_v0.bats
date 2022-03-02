#!/usr/bin/env bats

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'

    DIR="$( cd "$( dirname "$BATS_TEST_FILENAME" )" >/dev/null 2>&1 && pwd )"
    PATH="$DIR/../:$PATH"
    TEST_TEMP_DIR="$(mktemp -d)"
    source alignment_free_tool_workflow.sh
    cd $DIR/..
    jellyfish_count $TEST_TEMP_DIR
    jellyfish_dump $TEST_TEMP_DIR
    cat ${TEST_TEMP_DIR}/DI-1-1_S6.CharFreq
    cal_d2s $TEST_TEMP_DIR
    generate_matrix $TEST_TEMP_DIR
}

teardown() {
    rm -rf "$TEST_TEMP_DIR"
}

check_jf() {
    cmp data/jf/DI-1-1_S6.jf ${TEST_TEMP_DIR}/DI-1-1_S6.jf
    cmp data/jf/FI-2-21_S28.jf ${TEST_TEMP_DIR}/FI-2-21_S28.jf
    cmp data/jf/MI-1-19_S9.jf ${TEST_TEMP_DIR}/MI-1-19-S9.jf
    cmp data/jf/TAY_9_S28.jf ${TEST_TEMP_DIR}/TAY_9_S28.jf
}

check_charfreq() {
    cmp data/charfreq/DI-1-1_S6.CharFreq ${TEST_TEMP_DIR}/DI-1-1_S6.CharFreq
    cmp data/charfreq/FI-2-21_S28.CharFreq ${TEST_TEMP_DIR}/FI-2-21_S28.CharFreq
    cmp data/charfreq/MI-1-19_S9.CharFreq ${TEST_TEMP_DIR}/MI-1-19-S9.CharFreq
    cmp data/charfreq/TAY_9_S28.CharFreq ${TEST_TEMP_DIR}/TAY_9_S28.CharFreq
}

check_nkz() {
    cmp data/nkz/DI-1-1_S6.nkz ${TEST_TEMP_DIR}/DI-1-1_S6.nkz
    cmp data/nkz/FI-2-21_S28.nkz ${TEST_TEMP_DIR}/FI-2-21_S28.nkz
    cmp data/nkz/MI-1-19_S9.nkz ${TEST_TEMP_DIR}/MI-1-19-S9.nkz
    cmp data/nkz/TAY_9_S28.nkz ${TEST_TEMP_DIR}/TAY_9_S28.nkz
}

check_d2s() {
    cmp data/d2s/DI-1-1_S6-FI-2-21_S28.txt   ${TEST_TEMP_DIR}/DI-1-1_S6-FI-2-21_S28.txt
    cmp data/d2s/DI-1-1_S6-MI-1-19_S9.txt    ${TEST_TEMP_DIR}/DI-1-1_S6-MI-1-19_S9.txt
    cmp data/d2s/DI-1-1_S6-TAY_9_S28.txt     ${TEST_TEMP_DIR}/DI-1-1_S6-TAY_9_S28.txt
    cmp data/d2s/FI-2-21_S28-MI-1-19_S9.txt  ${TEST_TEMP_DIR}/FI-2-21_S28-MI-1-19_S9.txt
    cmp data/d2s/FI-2-21_S28-TAY_9_S28.txt   ${TEST_TEMP_DIR}/FI-2-21_S28-TAY_9_S28.txt
    cmp data/d2s/MI-1-19_S9-TAY_9_S28.txt    ${TEST_TEMP_DIR}/MI-1-19_S9-TAY_9_S28.txt
}

check_matrix() {
    cmp data/matrix.txt ${TEST_TEMP_DIR}/matrix.txt
}


@test "the output from jellyfish count" {
    run check_jf
    [ "$output" = "" ]
}

@test "charfreq output" {
    run check_charfreq
    [ "$output" = "" ] 
}

@test "nkz output" {
    run check_nkz
    [ "$output" = "" ] 
}

@test "the output from d2s" {
    run check_d2s
    [ "$output" = "" ]
}

@test "the final matrix" {
    run check_matrix
    [ "$output" = "" ]
}

