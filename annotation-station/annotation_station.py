import argparse
import os
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('input_file', type=str,
        help='input file')

# example of an argument group
input_str_group = parser.add_argument_group('input_str_group')
input_str_group.add_argument('--input-string', type=str,
        help='an input string')

parser.add_argument('input_bool', action='store_true',
        help='fp to be used for readcount output')
parser.add_argument('--input-int', type=int,
        default=2, help='an input int')

args = parser.parse_args()

def check_arguments():
    if args.input_str is None:
        raise ValueError('Must specify an input string')


def run_stuff(input_file, input_str, input_bool, input_int):
    """Run mytool"""
    tool_args = ['blah', 'beep',
            '-@', str(input_int),
            '-i', input_str]

    if input_bool:
        tool_args += ['-b', 'stuff']

    tool_args.append(input_file)

    print('mytool is executing blah step')
    print(f'blah is executing the following command: {" ".join(tool_args)}')
    print(subprocess.check_output(tool_args).decode('utf-8'))


def main():
    check_arguments()

    run_stuff(args.input_file, args.input_str, args.input_book, args.input_int)

if __name__ == '__main__':
    main()
