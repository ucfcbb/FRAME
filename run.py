import argparse
from inference import Chrom

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="input IBD matches per chrom")
    parser.add_argument("--r", required=True, help="reference population file")
    parser.add_argument("--q", required=True, help="query sample file")
    parser.add_argument("--o", required=True, help="output results folder")
    parser.add_argument("--n", required=True, help="numSitesPerChrom file")
    parser.add_argument("--c", required=True, help="target chromosome")

    args = parser.parse_args()
    input_ibd_file = args.i
    ref_popln_file = args.r
    query_sample_file = args.q
    output_folder = args.o
    numSiteFile = args.n
    chrom_num = args.c

    chromosome = Chrom.Chrom(ref_popln_file, query_sample_file, output_folder, 
                        numSiteFile, chrom_num)
    chromosome.process_all_matches(input_ibd_file)

if __name__ == '__main__':
    main()
