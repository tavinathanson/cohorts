"""
Utility functions for tests.
"""
from os import path

DATA_DIR = path.join(path.dirname(__file__), "data")

def data_path(name):
    """
    Return the absolute path to a file in the varcode/test/data directory.
    The name specified should be relative to varcode/test/data.
    """
    return path.join(DATA_DIR, name)

def generated_data_path(name):
    return path.join(DATA_DIR, "generated", name)
