import argparse
import os
from natsort import natsorted
from inference import Chrom
from inference import helper

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="input IBD matches dir ")
    parser.add_argument("--r", required=True, help="reference population file")
    parser.add_argument("--q", required=True, help="query sample file")
    parser.add_argument("--o", required=True, help="output results folder")
    parser.add_argument("--n", required=True, help="numSitesPerChrom file")

    args = parser.parse_args()
    input_ibd_dir = args.i
    ref_popln_file = args.r
    query_sample_file = args.q
    output_folder = args.o
    numSiteFile = args.n
    chrom_num = args.c

    # calculate weights per chrom
    ibd_files = os.listdir(input_ibd_dir)
    for num, fl in enumerate(natsorted(ibd_files)):
        chrom_num = num + 1
        input_ibd_file = os.path.join(input_ibd_dir, fl)
        assert os.path.exists(input_ibd_file)
        chromosome = Chrom.Chrom(ref_popln_file, query_sample_file, output_folder, 
                            numSiteFile, chrom_num)
        chromosome.process_all_matches(input_ibd_file)

    # process all weights
    helper.process_all_weights(numSiteFile, output_folder)

    # calculate proportions
    helper.calculate_proportions('all-chr-ancestry-normed.json')

if __name__ == '__main__':
    main()
