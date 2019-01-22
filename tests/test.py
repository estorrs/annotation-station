import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'
TEST_INPUT_FILE_1 = os.path.join(TEST_DATA_DIR, 'threshold.tsv')
TEST_INPUT_FILE_2 = os.path.join(TEST_DATA_DIR, 'test.tsv')
TEST_OUTPUT_FILE_1 = os.path.join(TEST_DATA_DIR, 'threshold.output.tsv')
TEST_OUTPUT_FILE_2 = os.path.join(TEST_DATA_DIR, 'test.output.tsv')

def test_transvar_annotation():
    tool_args = ['python', 'annotation-station/annotation_station.py',
            '--input-header',
            '--output', TEST_OUTPUT_FILE_1,
            '--input-type', 'tsv',
            TEST_INPUT_FILE_1]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert True

def test_transvar_annotation_brca():
    tool_args = ['python', 'annotation-station/annotation_station.py',
            '--input-header',
            '--output', TEST_OUTPUT_FILE_2,
            '--input-type', 'tsv',
            TEST_INPUT_FILE_2]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert True
