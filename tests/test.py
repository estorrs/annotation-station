import os
import subprocess

import pytest

TEST_DATA_DIR = 'tests/data/'
TEST_INPUT_FILE = os.path.join(TEST_DATA_DIR, 'test_file.txt')

def test_mytool():
    tool_args = ['python', 'mytool/mytool.py',
            '--input-int', '2',
            '--input_string', 'hello',
            '--input-bool',
            TEST_INPUT_FILE]
    
    results = subprocess.check_output(tool_args).decode('utf-8')

    assert True
