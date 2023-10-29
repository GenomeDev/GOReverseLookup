"""Script for preforming reverse lookup study with an input file"""

import argparse
from goreverselookup.Model import InputFileParser, ReverseLookupModel


def main():
    parser = argparse.ArgumentParser(description="Run GO reverse lookup study")
    parser.add_argument("filename", help="the filepath to the input file")
    parser.add_argument(
        "-v",
        "--verbose",
        help="change the verbosity of console output",
        action="count",
        default=0,
        type=int,
    )
    parser.add_argument("-o", "--output", help="destination filepath for the report")
    parser.add_argument("--noapi", help="disable the use of APIs in the study")

    args = parser.parse_args()

    input_file = InputFileParser(args.filename)

    model = ReverseLookupModel(
        goterms_per_process=
        target_processes=
        obo_filepath=
        annotations_filepaths=
        target_species=
        ortholog_species=
        valid_relationships=
        valid_evidence_codes=
        subontologies=
        indirect_annotations_propagation=
        alpha=
        pvalcal=
        correction_methods=
    )

    # after sucessfully reading the whole file, parse data from supplied files
    # TODO: move most of this calls to wrapers in methods, so it is more pythonic
    # TODO: move this part into the script file, it will give better control

    model.godag = GODag(model.obo_filepath)
    # TODO: filter subontologies
    model.godag.filter(lambda a: a.namespace in model.subontologies)

    model.annotations = Annotations().union(Annotations.from_file(fp) for fp in model.annotations_filepaths)
    # TODO: filter evidence codes

    # filter annotations to only leave the ones in ortholog_species
    model.annotations.filter(lambda a: a.taxon in [model.target_species, *model.ortholog_species])

    # find orthologs
    model.find_orthologs(target_species=model.target_species)
    # discard all annotations for which an ortholog was not found
    model.annotations.filter(lambda a: a.taxon == model.target_species)

    # convert to ensg
    model.annotations.convert_ids("ensg", "gConvert")

    # we don't need genes, which are not associated to any goterm
    goterm_list = [gt.id for gt in model.goterms]
    all_relevant_genes = [object_id for object_id, anno_set in model.annotations.dict_from_attr("object_id").items() if any(anno.term_id in goterm_list for anno in anno_set)]
    model.annotations.filter(lambda a: a.object_id in all_relevant_genes)

    # propagate associations
    # TODO different types
    if model.indirect_annotations_propagation is not False:
        model.annotations.propagate_associations(model.godag)



if __name__ == "__main__":
    main()
