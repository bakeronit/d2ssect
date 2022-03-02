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
    #cal_d2s $TEST_TEMP_DIR
    #generate_matrix $TEST_TEMP_DIR
}

teardown() {
    rm -rf "$TEST_TEMP_DIR"
}

check_charfreq() {
    cmp data/DI-1-1_S6.CharFreq ${TEST_TEMP_DIR}/DI-1-1_S6.CharFreq
}

#@test "the output from jellyfish should be the same, this is a binary compressed file" 


@test "the output from charfreq should be the same, text file" {
    run check_charfreq
    [ "$output" = "" ] 
}

#@test "the output from d2s, a single text file"

#@test "the final matrix, text file"

