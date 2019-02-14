import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'

TEST_GENE_TO_PRIMARY_TRANSCRIPT_FP = os.path.join(TEST_DATA_DIR,
        'test.gene_to_primary_transcript.tsv')
TEST_REPEATS_TABLE_FP = os.path.join(TEST_DATA_DIR,
        'test.repeats_table.tsv')

TEST_INPUT_FILE_1 = os.path.join(TEST_DATA_DIR, 'threshold.tsv')
TEST_INPUT_FILE_2 = os.path.join(TEST_DATA_DIR, 'test.tsv')
REPEATS_INPUT_FILE = os.path.join(TEST_DATA_DIR, 'repeats.tsv')
TEST_OUTPUT_FILE_1 = os.path.join(TEST_DATA_DIR, 'threshold.output.tsv')
TEST_OUTPUT_FILE_2 = os.path.join(TEST_DATA_DIR, 'test.output.tsv')
REPEATS_OUTPUT_FILE = os.path.join(TEST_DATA_DIR, 'repeats.output.tsv')

def test_transvar_annotation():
    tool_args = ['python', 'annotation-station/annotation_station.py',
            '--input-header',
            '--annotate-transvar',
            '--primary-transcripts', TEST_GENE_TO_PRIMARY_TRANSCRIPT_FP,
            '--output', TEST_OUTPUT_FILE_1,
            '--input-type', 'tsv',
            TEST_INPUT_FILE_1]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert 'WASH7P' in open(TEST_OUTPUT_FILE_1).read()

def test_transvar_annotation_brca():
    tool_args = ['python', 'annotation-station/annotation_station.py',
            '--input-header',
            '--annotate-transvar',
            '--primary-transcripts', TEST_GENE_TO_PRIMARY_TRANSCRIPT_FP,
            '--output', TEST_OUTPUT_FILE_2,
            '--input-type', 'tsv',
            TEST_INPUT_FILE_2]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    o = open(TEST_OUTPUT_FILE_2).read()
    assert 'BRCA1' in o and 'ENST00000471181' in o

def test_repeats_annotation():
    tool_args = ['python', 'annotation-station/annotation_station.py',
            '--input-header',
            '--annotate-repeats',
            '--repeats-table', TEST_REPEATS_TABLE_FP,
            '--output', REPEATS_OUTPUT_FILE,
            '--input-type', 'tsv',
            REPEATS_INPUT_FILE]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    l = [x for x in open(REPEATS_OUTPUT_FILE) if 'AluSc' in x][0]
    assert '43048295' in l

def test_all_annotation():
    tool_args = ['python', 'annotation-station/annotation_station.py',
            '--input-header',
            '--annotate-repeats',
            '--repeats-table', TEST_REPEATS_TABLE_FP,
            '--annotate-transvar',
            '--primary-transcripts', TEST_GENE_TO_PRIMARY_TRANSCRIPT_FP,
            '--output', REPEATS_OUTPUT_FILE,
            '--input-type', 'tsv',
            REPEATS_INPUT_FILE]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    l = [x for x in open(REPEATS_OUTPUT_FILE) if 'AluSc' in x][0]
    assert '43048295' in l and 'BRCA1' in l and 'ENST00000471181' in l
