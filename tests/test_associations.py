import pytest

from goreverselookup.associations import Annotation, Annotations
from goreverselookup.ontology import GODag, GOTerm


def test_header(annotations_test):
    assert annotations_test.version == "2.2"
    assert annotations_test.date == "2023-07-29T02:43"


@pytest.mark.skip
def test_UniProtKBA0A024RBG1_assoc(annotations_test):
    assert "GO:0002250" in annotations_test
    assert len(annotations_test["GO:0002250"]) == 2
    assert (
        next(
            obj.relationship
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "involved_in"
    )
    assert (
        next(
            obj.NOTrelation
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        is False
    )
    assert (
        next(
            obj.reference
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "GO_REF:0000043"
    )
    assert (
        next(
            obj.evidence_code
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "IEA"
    )
    assert (
        next(
            obj.taxon
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "taxon:9606"
    )
    assert (
        next(
            obj.date
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "20230703"
    )


@pytest.mark.skip
def test_UniProtKBA0A024RBG1_cardinality_0_fields_assoc(annotations_test):
    assert (
        next(
            obj.relationship
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "located_in"
    )
    assert (
        next(
            obj.NOTrelation
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        is False
    )
    assert (
        next(
            obj.reference
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "GO_REF:0000052"
    )
    assert (
        next(
            obj.evidence_code
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "IDA"
    )
    assert (
        next(
            obj.taxon
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "taxon:9606"
    )
    assert (
        next(
            obj.date
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "20230619"
    )


def test_propagate_associations(annotations_test: Annotations, godag_test: GODag):
    annotations_test.propagate_associations(godag_test)
    assert any(anno.term_id == "GO:0000001" for anno in annotations_test)
    assert not any(anno.term_id == "GO:0000003" for anno in annotations_test)
    assert sum(1 for anno in annotations_test if anno.term_id == "GO:0000001") == 2


def test_dict_from_attr():
    anno1 = Annotation(object_id="ABC1", term_id="GO:1234")
    anno2 = Annotation(object_id="ABC2", term_id="GO:1234")
    anno3 = Annotation(object_id="ABC2", term_id="GO:5678")
    annoset = Annotations([anno1, anno2, anno3])

    dict_by_term_id = annoset.dict_from_attr("term_id")
    assert len(dict_by_term_id) == 2
    assert len(dict_by_term_id["GO:1234"]) == 2

    dict_by_object_id = annoset.dict_from_attr("object_id")
    assert len(dict_by_object_id) == 2
    assert len(dict_by_object_id["ABC2"]) == 2


def test_annotation_set_operations():
    anno1 = Annotation(object_id="ABC1", term_id="GO:1234")
    anno2 = Annotation(object_id="ABC2", term_id="GO:1234")
    anno3 = Annotation(object_id="ABC2", term_id="GO:5678")

    annoset1 = Annotations([anno1])
    annoset2 = Annotations([anno2])
    annoset3 = Annotations([anno3])

    assert annoset1.union(annoset2) == annoset2.union(
        annoset1
    )  # union should be the same regardles of order

    assert isinstance(annoset1.union(annoset2), Annotations)

    assert len(annoset1.union(annoset2, annoset3)) == 3

    assert annoset1.intersection(annoset2) == Annotations()

    assert annoset1.intersection(annoset1.union(annoset2)) == annoset1

    with pytest.raises(TypeError):
        annoset1.union({"a", "b"})


def test_match_annotations_to_godag():
    annoset = Annotations(
        [
            Annotation(object_id="ABC1", term_id="GO:1234"),
            Annotation(object_id="ABC2", term_id="GO:1234"),
            Annotation(object_id="ABC2", term_id="GO:5678"),
        ]
    )

    godag = GODag()
    godag["GO:1234"] = GOTerm("GO:1234")

    assert len(annoset) == 3

    annoset.match_annotations_to_godag(godag)
    assert len(annoset) == 2


def test_add_taxon_to_object_id():
    annoset = Annotations(
        [
            Annotation(object_id="ABC1", term_id="GO:1234", taxon="9606"),
            Annotation(object_id="ABC2", term_id="GO:1234"),
        ]
    )
    annoset.add_taxon_to_object_id()

    assert any(a.object_id == "ABC1-9606" for a in annoset)
    assert any(a.object_id == "ABC2" for a in annoset)


def test_find_orthologs_gOrth():
    anno1 = Annotation(
        object_id="ZFIN:ZDB-GENE-040912-6", term_id="GO:1234", taxon="7955"
    )  # returns multiple
    anno2 = Annotation(
        object_id="ZFIN:ZDB-GENE-170217-1", term_id="GO:1234", taxon="7955"
    )  # returns one
    anno2_n = Annotation(
        object_id="ZFIN:ZDB-GENE-170217-1", term_id="GO:5678", taxon="7955"
    )  # returns one
    anno3 = Annotation(
        object_id="ZFIN:ZDB-GENE-021119-1", term_id="GO:5678", taxon="7955"
    )  # returns two, but one N/A
    annoset = Annotations([anno1, anno2, anno2_n, anno3])

    annoset.find_orthologs(taxon="9606", database="gOrth", prune=True)

    assert len(annoset) == 6
    assert all(r.taxon == "9606" for r in annoset)  # test if prune works
    assert set(a.object_id for a in annoset) == set(
        [
            "ENSG00000005421",
            "ENSG00000105852",
            "ENSG00000105854",
            "ENSG00000151577",
            "ENSG00000168938",
        ]
    )


def test_filter():
    annoset = Annotations(
        [
            Annotation(object_id="ABC1", term_id="GO:1234", taxon="9606"),
            Annotation(object_id="ABC2", term_id="GO:1234"),
        ]
    )
    annoset.filter(lambda a: a.taxon == "9606")
    assert len(annoset) == 1
    assert Annotation(object_id="ABC1", term_id="GO:1234", taxon="9606") in annoset


def test_convert_ids_gConvert():
    anno1 = Annotation(
        object_id="ZFIN:ZDB-GENE-040912-6", term_id="GO:1234", taxon="7955"
    )  # returns multiple
    anno2 = Annotation(
        object_id="ZFIN:ZDB-GENE-170217-1", term_id="GO:1234", taxon="7955"
    )  # returns one
    anno2_n = Annotation(
        object_id="ZFIN:ZDB-GENE-170217-1", term_id="GO:5678", taxon="7955"
    )  # returns one
    anno3 = Annotation(
        object_id="ZFIN:ZDB-GENE-021119-1", term_id="GO:5678", taxon="7955"
    )  # returns two, but one N/A
    annoset = Annotations([anno1, anno2, anno2_n, anno3])

    annoset.convert_ids(namespace="ensg", database="gConvert")

    assert len(annoset) == 5
    assert set(a.object_id for a in annoset) == set(
        [
            "ENSDARG00000032496",
            "ENSDARG00000032131",
            "ENSDARG00000110679",
            "ZFIN:ZDB-GENE-170217-1",
        ]
    )
