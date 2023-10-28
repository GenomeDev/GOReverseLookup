"""Script for preforming reverse lookup study with an input file"""

import argparse

def main():
    parser = argparse.ArgumentParser(description='Run GO reverse lookup study')
    parser.add_argument('input', help="the filepath to the input file")
    parser.add_argument('-v', '--verbose', help="change the verbosity of console output", action='count', default=0, type=int)
    parser.add_argument('-o', '--output', help="destination filepath for the report")
    parser.add_argument('--noapi', help="disable the use of APIs in the study")

    args = parser.parse_args()

    

if __name__ == "__main__":
    main()