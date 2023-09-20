import argparse
import os
import re
from natsort import natsorted
from inference import Chrom
from inference import helper

# pybind11
import sys

#sys.path.append("build/")
import ibd_call


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=
        "input dir with (reference and query VCF files by default, unless used with --precomputed"
    )

    parser.add_argument("-l",
                        "--length",
                        required=True,
                        help="length of IBD matches (in sites)")

    parser.add_argument("-p",
                        "--refpop",
                        required=True,
                        help="reference population file")
    parser.add_argument("-s",
                        "--sample",
                        required=True,
                        help="query sample file")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="output results folder")
    parser.add_argument("-n",
                        "--sitesperchrom",
                        required=True,
                        help="numSitesPerChrom file")

    # optional arguments
    parser.add_argument("-B",
                        "--bits",
                        type=int,
                        choices=[64, 128],
                        default=128,
                        help="bits size: 64 or 128 (default)")
    parser.add_argument(
        "-c",
        "--chrom",
        type=int,
        choices=list(range(1, 23)),
        help="chromosome number (individual chromosome to be processed)")
    parser.add_argument("--precomputed",
                        help="pre-computed input IBD matches",
                        action="store_true")

    args = parser.parse_args()

    match_length = int(args.length)
    ref_popln_file = args.refpop
    query_sample_file = args.sample
    output_folder = args.output
    numSiteFile = args.sitesperchrom

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # if individual chromosome to be processed
    targetChrom = None
    if args.chrom:
        targetChrom = args.chrom

    # use pre-computed IBDs, if they exist
    if args.precomputed:
        # add sanity check to see that matches exist
        print("Using precomputed IBDs...", flush=True)
        input_ibd_dir = args.input
    else:
        # call IBDs using Syllable-PBWT
        input_vcfs_dir = args.input
        if (args.bits):
            bits = int(args.bits)
        else:
            bits = 128

        query_vcfs = [
            f for f in natsorted(os.listdir(input_vcfs_dir)) if 'query' in f
        ]
        ref_vcfs = [
            f for f in natsorted(os.listdir(input_vcfs_dir)) if 'ref' in f
        ]

        for q_vcf, r_vcf in zip(query_vcfs, ref_vcfs):
            query_file = os.path.join(input_vcfs_dir, q_vcf)
            ref_file = os.path.join(input_vcfs_dir, r_vcf)
            temp = re.findall(r'\d+', q_vcf)
            chrom = list(map(int, temp))[0]

            if targetChrom and chrom == targetChrom:
                print(f"Calculating IBDS for only CHR-{chrom}", flush=True)
                ibd_call.call_ibds(ref_file, query_file, output_folder,
                                   match_length, bits, str(chrom))
            else:
                print(f"Calculating IBDS for CHR-{chrom}", flush=True)

                ibd_call.call_ibds(ref_file, query_file, output_folder,
                                   match_length, bits, str(chrom))

        input_ibd_dir = os.path.join(output_folder, "sq_ibds")

    # calculate weights per chrom
    ibd_files = os.listdir(input_ibd_dir)
    for num, fl in enumerate(natsorted(ibd_files)):
        #chrom_num = num + 1
        temp = re.findall(r'\d+', fl)
        chrom_num = list(map(int, temp))[0]
        input_ibd_file = os.path.join(input_ibd_dir, fl)
        print(
            f"Chr {chrom_num}: Calculating ancestry dosages: {input_ibd_file} ",
            flush=True)
        assert os.path.exists(input_ibd_file)
        chromosome = Chrom.Chrom(ref_popln_file, query_sample_file,
                                 output_folder, numSiteFile, str(chrom_num))
        chromosome.process_all_matches(input_ibd_file)

    # process all weights
    weights_dir = os.path.join(output_folder, 'weights')
    helper.process_all_weights(numSiteFile, weights_dir)

    # calculate proportions
    helper.calculate_proportions(
        os.path.join(weights_dir, 'all-chr-ancestry-normed.json'),
        output_folder)


if __name__ == '__main__':
    main()
