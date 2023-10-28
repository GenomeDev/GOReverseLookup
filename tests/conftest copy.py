import os

import pytest

from revonto.associations import Annotations
from revonto.ontology import GODag


@pytest.fixture
def annotations_test():
    return Annotations.from_file(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf")
    )


@pytest.fixture
def godag_test():
    return GODag.from_file(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/go1.obo")
    )
