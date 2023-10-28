from goreverselookup.ontology import GODag


def test_obo_dataversion(godag_test: GODag):
    assert godag_test.data_version == "test/2023-09-03"


def test_number_of_elements(godag_test: GODag):
    assert len(godag_test) == 6  # 6 non-obsolete terms + 1


def test_GO0000006_entry(godag_test: GODag):
    entry = godag_test["GO:0000006"]
    assert entry.term_id == "GO:0000006"
    assert entry.name == "third level 2"
    assert entry.namespace == "molecular_function"
    assert entry.description == '"child of second level 1" [TC:2.A.5.1.1]'
    assert entry._parents == {"GO:0000002"}
    assert entry.parents == {godag_test["GO:0000002"]}  # direct parents
    assert entry.children == {godag_test["GO:0000015"]}  # direct children
    assert entry.depth == 2
    assert entry.height == 1


def test_GO0000015_partof(godag_test: GODag):
    entry = godag_test["GO:0000015"]
    assert entry.parents == {godag_test["GO:0000006"], godag_test["GO:0005829"]}


def test_has_parent(godag_test: GODag):
    entry = godag_test["GO:0000002"]
    assert entry.has_parent("GO:0000001")


def test_has_child(godag_test: GODag):
    entry = godag_test["GO:0000002"]
    assert entry.has_child("GO:0000006")


def test_get_all_parents(godag_test: GODag):
    entry = godag_test["GO:0000006"]
    assert entry.get_all_parents() == {"GO:0000002", "GO:0000001"}


def test_get_all_children(godag_test: GODag):
    entry = godag_test["GO:0000002"]
    assert entry.get_all_children() == {"GO:0000006", "GO:0000015"}
