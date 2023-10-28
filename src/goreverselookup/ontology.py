"""
Read and store Gene Ontology's obo file.
Part of code has been taken from H Tang et al. 2018 (https://github.com/tanghaibao/goatools)
"""
# -*- coding: UTF-8 -*-
import os
from typing import Optional, Set

# if TYPE_CHECKING:
#    from .Metrics import Metrics, basic_mirna_score


class OBOReader(object):
    """Read goatools.org's obo file. Load into this iterable class.

    Download obo from: http://geneontology.org/ontology/go-basic.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print(rec)
    """

    # Scalar attributes for Typedefs:
    #                    'is_class_height', 'is_metadata_tag',
    #                    'is_transitive', 'transitive_over'])

    def __init__(self, obo_file="go-basic.obo"):
        """Read obo file. Load dictionary."""
        self.format_version = None  # e.g., "1.2" of "format-version:" line
        self.data_version = (
            None  # e.g., "releases/2016-07-07" from "data-version:" line
        )

        # True if obo file exists or if a link to an obo file exists.
        if os.path.isfile(obo_file):
            self.obo_file = obo_file
            # GOTerm attributes that are necessary for any operations:
        else:
            raise Exception(
                "COULD NOT READ({OBO})\n"
                "download obo file first\n "
                "[http://geneontology.org/ontology/"
                "go-basic.obo]".format(OBO=obo_file)
            )

    def __iter__(self):
        """Return one GO Term record at a time from an obo file."""
        # Wait to open file until needed. Automatically close file when done.
        with open(self.obo_file) as fstream:
            hdr = True
            rec_curr = None  # Stores current GO Term
            typedef_curr = None  # Stores current typedef
            for line in fstream:
                # obo lines start with any of: [Term], [Typedef], /^\S+:/, or /^\s*/
                if hdr:
                    if not self._init_obo_hdr(line):
                        hdr = False
                if rec_curr is None and line[0:6].lower() == "[term]":
                    rec_curr = GOTerm()
                elif typedef_curr is None and line[0:9].lower() == "[typedef]":
                    typedef_curr = True
                elif rec_curr is not None or typedef_curr is not None:
                    line = line.rstrip()  # chomp
                    if line:
                        self._add_to_obj(rec_curr, typedef_curr, line)
                    else:
                        if rec_curr is not None:
                            yield rec_curr
                            rec_curr = None
                        elif typedef_curr is not None:
                            typedef_curr = None
            # Return last record, if necessary
            if rec_curr is not None:
                yield rec_curr

    def _add_to_obj(self, rec_curr, typedef_curr, line):
        """Add information on line to GOTerm or Typedef."""
        if rec_curr is not None:
            self._add_to_ref(rec_curr, line)
        else:
            pass

    def _init_obo_hdr(self, line):
        """Save obo version and release."""
        if line[0:14] == "format-version":
            self.format_version = line[16:-1]
            return True
        if line[0:12] == "data-version":
            self.data_version = line[14:-1]
            return True
        if line[:17] == "default-namespace":
            self.default_namespace = line[18:].strip()
            return True
        if line[0:6].lower() == "[term]":
            return False
        return True

    def _add_to_ref(self, rec_curr, line):
        """Add new fields to the current reference."""
        # Examples of record lines containing ':' include:
        #   id: GO:0000002
        #   name: mitochondrial genome maintenance
        #   namespace: biological_process
        #   def: "The maintenance of ...
        #   is_a: GO:0007005 ! mitochondrion organization
        if line[:4] == "id: ":
            assert not rec_curr.term_id
            term_id = line[4:]
            rec_curr.term_id = term_id
            rec_curr.id = term_id
        elif line[:8] == "alt_id: ":
            rec_curr.alt_ids.add(line[8:])
        elif line[:6] == "name: ":
            assert not rec_curr.name
            rec_curr.name = line[6:]
        elif line[:5] == "def: ":
            rec_curr.description = line[5:]
        elif line[:11] == "namespace: ":
            rec_curr.namespace = line[11:]
        elif (
            line[:6] == "is_a: "
        ):  # based on https://geneontology.org/docs/ontology-relations/ "is_a" or "part_of" can be safely used as group annotations
            rec_curr._parents.add(line[6:].split()[0])
        elif line[:22] == "relationship: part_of ":
            rec_curr._parents.add(line[22:].split()[0])
        elif line[:13] == "is_obsolete: " and line[13:] == "true":
            rec_curr.is_obsolete = True


class GOTerm(object):
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(
        self,
        term_id: str = "",
        name: str = "",
        description: str = "",
        namespace: str = "default",
        is_obsolete: bool = False,
    ):
        self.term_id = term_id  # GO:NNNNNNN
        self.name = name  # description
        self.description = description
        self.namespace = namespace  # BP, CC, MF
        self._parents: set[str] = set()  # is_a basestring of parents
        self.parents: set[GOTerm] = set()  # direct parent records
        self.children: set[GOTerm] = set()  # direct children records
        self.is_obsolete = is_obsolete  # is_obsolete
        self.alt_ids: set[str] = set()  # alternative identifiers
        self.height: Optional[int] = None
        self.depth: Optional[int] = None

    def has_parent(self, term):
        """Return True if this GO object has a parent GO ID."""
        for parent in self.parents:
            if parent.term_id == term or parent.has_parent(term):
                return True
        return False

    def has_child(self, term):
        """Return True if this GO object has a child GO ID."""
        for parent in self.children:
            if parent.term_id == term or parent.has_child(term):
                return True
        return False

    def get_all_parents(self):
        """Return all parent GO IDs."""
        all_parents = set()
        for parent in self.parents:
            all_parents.add(parent.term_id)
            all_parents |= parent.get_all_parents()
        return all_parents

    def get_all_children(self) -> Set[str]:
        """Return all child GO IDs."""
        all_children = set()
        for child in self.children:
            all_children.add(child.term_id)
            all_children |= child.get_all_children()
        return all_children


class GODag(dict[str, GOTerm]):
    """Holds the GO DAG as a dict."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_file(cls, file, load_obsolete=False):
        """Read obo file. Store results."""
        reader = OBOReader(file)

        instance = cls()
        # TODO: Save alt_ids and their corresponding main GO ID. Add to GODag after populating GO Terms
        for rec in reader:
            # Save record if:
            #   1) Argument load_obsolete is True OR
            #   2) Argument load_obsolete is False and the GO term is "live" (not obsolete)
            if load_obsolete or not rec.is_obsolete:
                instance[rec.term_id] = rec

        instance._populate_terms()
        instance._set_height_depth()

        # TODO: Add alt_ids to go2obj
        # for goid_alt, rec in alt2rec.items():
        #    self[goid_alt] = rec
        desc = instance._str_desc(reader)

        instance.version = desc
        instance.data_version = reader.data_version

        return instance

    def _str_desc(self, reader):
        """String containing information about the current GO DAG."""
        data_version = reader.data_version
        if data_version is not None:
            data_version = data_version.replace("releases/", "")
        desc = "{OBO}: fmt({FMT}) rel({REL}) {N:,} Terms".format(
            OBO=reader.obo_file,
            FMT=reader.format_version,
            REL=data_version,
            N=len(self),
        )
        return desc

    def _populate_terms(self):
        """Convert GO IDs to GO Term record objects. Populate children."""

        # Make parents and relationships references to the actual GO terms.
        for rec in self.values():
            # Given parent GO IDs, set parent GO Term objects
            rec.parents = set(self[goid] for goid in rec._parents)

            # For each parent GO Term object, add it's child GO Term to the children data member
            for parent_rec in rec.parents:
                parent_rec.children.add(rec)

    def _set_height_depth(self):
        """Set height, depth and add inverted relationships."""

        def _init_height(rec: GOTerm) -> int:
            if rec.height is None:
                if rec.children:
                    rec.height = max(_init_height(rec) for rec in rec.children) + 1
                else:
                    rec.height = 0
            return rec.height

        def _init_depth(rec: GOTerm) -> int:
            if rec.depth is None:
                if rec.parents:
                    rec.depth = max(_init_depth(rec) for rec in rec.parents) + 1
                else:
                    rec.depth = 0
            return rec.depth

        for rec in self.values():
            if rec.height is None:
                _init_height(rec)

            if rec.depth is None:
                _init_depth(rec)

    @staticmethod
    def id2int(go_id):
        """Given a GO ID, return the int value."""
        return int(go_id.replace("GO:", "", 1))
