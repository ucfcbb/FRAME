import argparse
import os
from natsort import natsorted
from inference import Chrom
from inference import helper

# pybind11
import sys
sys.path.append("build/")
import ibd_call

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input", help="input dir with (reference and query VCF files by default, unless used with --precomputed")
    parser.add_argument("-L","--length", required=True, help="length of IBD matches (in sites)")
    parser.add_argument("-p","--refpop", required=True, help="reference population file")
    parser.add_argument("-s","--sample", required=True, help="query sample file")
    parser.add_argument("-o","--output", required=True, help="output results folder")
    parser.add_argument("-n","--sitesperchrom", required=True, help="numSitesPerChrom file")

    # optional arguments
    parser.add_argument("-B","--bits", type=int, choices=[64, 128], default=128, help="bits size: 64 or 128 (default)")
    parser.add_argument("--precomputed", help="pre-computed input IBD matches", action="store_true")

    args = parser.parse_args()

    match_length = int(args.L)
    ref_popln_file = args.p
    query_sample_file = args.s
    output_folder = args.o
    numSiteFile = args.n

    if (args.precomputed):
        # add sanity check to see that matches exist
        input_ibd_dir = args.i
    else:
        # call IBDs using Syllable-PBWT
        input_vcfs_dir = args.i
        if (args.B):
            bits = int(args.B)
        else:
            bits = 128
        
        query_vcfs = [f for f in os.listdir(input_vcfs_dir) if 'query' in f]
        ref_vcfs = [f for f in os.listdir(input_vcfs_dir) if 'ref' in f]

        for q_vcf, r_vcf in natsorted(zip(query_vcfs, ref_vcfs)):
            query_file = os.path.join(input_vcfs_dir, q_vcf)
            ref_file = os.path.join(input_vcfs_dir, r_vcf)
            ibd_call.call_ibds(ref_file, query_file, match_length, bits)


    # calculate weights per chrom
    ibd_files = os.listdir(input_ibd_dir)
    for num, fl in enumerate(natsorted(ibd_files)):
        chrom_num = num + 1
        input_ibd_file = os.path.join(input_ibd_dir, fl)
        assert os.path.exists(input_ibd_file)
        chromosome = Chrom.Chrom(ref_popln_file, query_sample_file, output_folder, 
                            numSiteFile, str(chrom_num))
        chromosome.process_all_matches(input_ibd_file)

    # process all weights
    helper.process_all_weights(numSiteFile, output_folder)

    # calculate proportions
    helper.calculate_proportions('all-chr-ancestry-normed.json')

if __name__ == '__main__':
    main()
