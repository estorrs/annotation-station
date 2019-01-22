import subprocess

def setup_transvar(reference_fasta_fp=None):
    """set up transvar database"""
    if reference_fasta_fp is None:
        tool_args = [

