import argparse
#import inference.

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="input batch folder")
    parser.add_argument("--rm", required=True, help="reference meta file")
    parser.add_argument("--r", required=True, help="reference vcf file")
    parser.add_argument("--o", required=True, help="output results folder")
    parser.add_argument("--n", required=True, help="numSitesPerChrom file")
    parser.add_argument("--b",
                        required=True,
                        help="batch processed flag (T/F)")
    args = parser.parse_args()
    input_folder = args.i
    output_folder = args.o
    ref_meta_file = args.rm
    ref_file = args.r
    numSiteFile = args.n
    batch_processed = args.b

    if batch_processed == 'F':
        chromosome = Chrom(input_folder, output_folder, ref_meta_file,
                           ref_file, numSiteFile, False)
        # input_folder is chromosome leve dir
        chromosome.process_all_matches(input_folder)
    else:
        chromosome = Chrom(input_folder, output_folder, ref_meta_file,
                           ref_file, numSiteFile, True)
        # input_folder is batch level dir
        chromosome.process_per_batch(input_folder)


if __name__ == '__main__':
    main()