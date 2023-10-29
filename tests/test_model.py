import pytest

from goreverselookup.Model import InputFileParser, ReverseLookupModel
import os

def test_InputFileParser():
    inputfile = InputFileParser(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/input.txt"))
    assert inputfile.obo_filepath == "data/go1.obo"
    assert inputfile.alpha == 0.01
    assert inputfile.subontologies == ["biological_process", "cellular_component", "molecular_function"]
    assert inputfile.ortholog_species == ["7955", "10116", "10090", "8353"]
    
def test_ReverseLookupModel_custom_data():
    model = ReverseLookupModel(
        obo_filepath=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/go1.obo"),
        annotations_filepaths=[os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf"), os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf")],
        subontologies=["biological_process", "cellular_component", "molecular_function"],
        goterms_per_process={"process1":{"GO:0000002", "GO:0005829"}})
    
    model.run_study()
    
    assert len(model.godag) == 6
    assert len(model.annotations) == 8
    assert model.results_dict["process1"][0].object_id == "UniProtKB:A0A024RBG1"
    assert pytest.approx(model.results_dict["process1"][0].pvals["uncorrected"]) == 0.20000000000000004
    assert (
        pytest.approx(model.results_dict["process1"][0].pvals["bonferroni"]) == 0.20000000000000004
    )
    
def test_ReverseLookupModel_from_input_file():
    inputfile = InputFileParser(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/input.txt"))
    inputfile.obo_filepath=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/go1.obo")
    inputfile.annotations_filepaths=[os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf")]
    model = ReverseLookupModel.from_input_file_parser(inputfile)
    
    model.run_study()
    
    assert len(model.godag) == 6
    assert len(model.annotations) == 8
    assert model.results_dict["chronic_inflammation+"][0].object_id == "UniProtKB:A0A024RBG1"
    assert pytest.approx(model.results_dict["chronic_inflammation+"][0].pvals["uncorrected"]) == 0.20000000000000004
    assert (
        pytest.approx(model.results_dict["chronic_inflammation+"][0].pvals["bonferroni"]) == 0.20000000000000004
    )