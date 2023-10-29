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
    parser.add_argument("--noapi", help="disable the use of APIs in the study", action="store_true")

    args = parser.parse_args()

    input_file = InputFileParser(args.filename)

    model = ReverseLookupModel.from_input_file_parser(input_file)

    # find orthologs
    if args.noapi:
        model.find_orthologs(database="local_files")
    else:
        model.find_orthologs(database="all")
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
        
    model.run_study()
    
    print(model.results_dict)

if __name__ == "__main__":
    main()
