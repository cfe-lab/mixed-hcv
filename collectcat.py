import sys
import os
import argparse
import subprocess
import re
from glob import glob
from datetime import datetime
from helper import Cache

def main():
    parser = argparse.ArgumentParser(
        description='Concatenates all runs together')
    
    parser.add_argument('-output', help='folder to output to', default=".")
    parser.add_argument('-runname', help='<input> name of the run')
    parser.add_argument('-x', help='path to bowtie2 index (.bt2)')
    parser.add_argument('-minq', type=int, help='minimum mapping quality (MAPQ)')
    parser.add_argument('-runtype', help='deli or full')
    parser.add_argument('-min_target_width', type=float, help='minimum fraction of target width covered to be considered a hit', default=0.0)
    parser.add_argument("--cache", help="the cache folder that holds all results and sam files", default=None)

    parser.add_argument('inputs', help='inputs', nargs='*', default=[])

    args = parser.parse_args()
    
    min_align_quality = args.minq
    min_target_width = args.min_target_width
    deli = args.runtype
    runname = args.runname
    cache_path = args.cache
    ref = os.path.basename(args.x)
    output = os.path.abspath(args.output)
    input_files = args.inputs

    # If we were given a cache
    if cache_path is not None:
        print "Found cache"
        # Try to dig up the results from there first
        cache = Cache(runname, min_align_quality, ref, deli, min_target_width, cache_path)
        files  = cache.list_cached_results()
        input_files = files


    # Make sure the filename is consistent
    if not os.path.isdir(os.path.join(output, runname)):
        os.makedirs(os.path.join(output, runname))

    filename = os.path.join(
        output, 
        runname,
        Cache._result_name(runname, deli, ref, min_align_quality, min_target_width)
    )

    # Concatenate all the inputs together
    with open(filename, "w") as target_file:
        header = ""
        results = []

        for csv in input_files:
            with open(csv, "r") as input_file:
                header = input_file.readline()
                results += input_file.readlines()

        target_file.write(header)
        target_file.write(''.join(results))


if __name__ == '__main__':
    main()