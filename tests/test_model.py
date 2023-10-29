import pytest

from goreverselookup.Model import InputFileParser
import os

def test_InputFileParser():
    inputfile = InputFileParser(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/input.txt"))
    assert inputfile.obo_filepath == "data/go1.obo"
    assert inputfile.alpha == 0.01
    assert inputfile.subontologies == ["BP", "MF"]
    assert inputfile.ortholog_species == ["7955", "10116", "10090", "8353"]