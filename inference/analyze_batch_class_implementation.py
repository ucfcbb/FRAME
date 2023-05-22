import os
import collections
import numpy as np
from natsort import natsorted
import argparse
import json
import re


class Chrom:
    # store all the calculated ancestral scores/weights
    __scores_per_site = None

    ##################################
    # Reference metadata information #
    ##################################
    # store ref-ancestry of each reference sample
    __ref_ancestry_of_sample_d = collections.defaultdict()

    # key: ref_ancestry, value: idx value
    __idx_of_ref_ancestry_d = collections.OrderedDict()
    __samples_per_ref_ancestry_d = collections.defaultdict(lambda: 0)

    # All the query sample information
    __num_sites_per_chrom = collections.defaultdict()

    def __init__(self,
                 input_dir,
                 output_folder,
                 ref_meta_file,
                 ref_file,
                 numSiteFile,
                 query_indv_file='indivs-list',
                 batch_processed=False):

        if not batch_processed:
            self.chrom_num = input_dir.split('/')[-1]
            print(
                f"-------------| Initializing {self.chrom_num} |-------------",
                flush=True)

            # setup output directory and output weights file
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)
            self.interm_output_path = os.path.join(
                output_folder,
                self.chrom_num + '-interm-ancestry-result-faster.json')
        else:
            batch = input_dir.split('/')[-1]
            self.batch_num = re.findall(r'\d+', batch)[0]
            self.chrom_num = input_dir.split('/')[-2]
            print(
                f"-------------| Initializing {self.chrom_num} |-------------",
                flush=True)
            print(f"    ---------|   Batch {self.batch_num}      |---------",
                  flush=True)

            # setup output directory and output weights file
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)
            intermediate_output_file = batch + '-interm-ancestry-result.json'
            self.interm_output_path = os.path.join(output_folder,
                                                   intermediate_output_file)

        # Load haplotype to Individual mapping for query & ref. samples
        self.haps2indivs_ref, self.uniq_indivs = self.get_indivs_ref(ref_file)
        self.haps2indivs_query = self.get_query_indivs(query_indv_file)

        # loading meta information
        self.get_num_sites_per_chrom(numSiteFile)
        self.load_all_required_meta_info(ref_meta_file, delim=" ")

        # Extra info
        self.total_num_sites = Chrom.__num_sites_per_chrom[self.chrom_num]
        self.total_num_uniq_ref_ancestries = len(Chrom.__idx_of_ref_ancestry_d)
        self.ref_ancestries = list(Chrom.__idx_of_ref_ancestry_d.keys())

        # create 2D numpy array : #sitesInChrom X #UniqRefAncestries
        Chrom.__scores_per_site = np.zeros(
            shape=(self.total_num_sites, self.total_num_uniq_ref_ancestries))

    def get_num_sites_per_chrom(self, numSiteFile):
        ''' 
            Stores the # of sites present in every chromosome's VCF file in a dict.
        '''
        with open(numSiteFile, 'r') as f:
            for line in f:
                if 'Chr' in line:
                    continue
                toks = line.strip().split('\t')
                Chrom.__num_sites_per_chrom[toks[0]] = int(toks[1])

    def get_query_hap_idx(self, input_file):
        '''
            returns the query haplotype's index
        '''
        return input_file.split('_')[-1]

    def get_query_indivs(self, query_indv_file):
        '''
            returns the mapping of query haplotype index to query individual id
        '''
        hap_indices = {}
        hap_id = 0
        #prefix = 'indivs-list'
        with open(query_indv_file, 'r') as f:
            for line in f:
                toks = line.rstrip().split()
                hap_indices[str(hap_id)] = toks[0]
                hap_id += 1
        return hap_indices

    def get_indivs_ref(self, ref_file):
        '''
            Make dictionary:
             key=ref-hap-idx, value= sample + '-0/1'
        '''
        haps2indivs_ref = {}
        uniq_indivs = set()
        with open(ref_file, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    toks = line.rstrip().split()[9:]
                    hap_id = 1
                    for t in toks:
                        haps2indivs_ref[str(hap_id)] = t + '-0'
                        haps2indivs_ref[str(hap_id + 1)] = t + '-1'
                        hap_id += 2
                        uniq_indivs.add(t)
                    break
        return haps2indivs_ref, uniq_indivs

    def load_all_required_meta_info(self, ref_meta_file, delim="\t"):
        idx = 0
        with open(ref_meta_file, 'r') as f:
            for line in f:
                if line.startswith('Sample'):
                    continue
                toks = line.rstrip().split(delim)
                sample_id = toks[0]
                population = toks[1]

                if sample_id not in self.uniq_indivs:
                    continue

                Chrom.__ref_ancestry_of_sample_d[sample_id] = population

                # store # of samples per ref_popln
                Chrom.__samples_per_ref_ancestry_d[population] += 1

                # assign some arbritrary index for each ref-popln
                # to initilize and index the scores_per_site array
                if population not in Chrom.__idx_of_ref_ancestry_d:
                    Chrom.__idx_of_ref_ancestry_d[population] = idx
                    idx += 1

    def process_all_matches(self, match_dir):
        ''' 
            Process all the matching segments for a chromosome ( w/o batching)
        '''
        match_files = os.listdir(match_dir)

        # create the interm output file
        fw_interm = open(self.interm_output_path, 'w')

        # iterate through all match files
        interim_res = {}
        normalized_weights = {}
        for ff in natsorted(match_files):
            print("-->{}".format(ff), flush=True)
            input_file = os.path.join(match_dir, ff)
            query_hap_idx = self.get_query_hap_idx(input_file)
            query_hap_id = self.haps2indivs_query[query_hap_idx]

            # Update #hits per site
            with open(input_file, 'r') as f:
                for line in f:
                    toks = line.rstrip().split()
                    ref_hap_idx = toks[2]
                    start_site = int(toks[3])
                    end_site = int(toks[4])  # exclusive

                    hit_ref_hap_idx = self.haps2indivs_ref[ref_hap_idx]
                    hit_ref_hap_sample = hit_ref_hap_idx[:-2]

                    array_idx_ref_ancestry = Chrom.__idx_of_ref_ancestry_d[
                        Chrom.__ref_ancestry_of_sample_d[hit_ref_hap_sample]]
                    Chrom.__scores_per_site[start_site:end_site,
                                            array_idx_ref_ancestry] += 1

            # process weights
            total_per_site = np.sum(Chrom.__scores_per_site, axis=1)
            for j in range(self.total_num_uniq_ref_ancestries):
                ref_anc_scores_across_sites = Chrom.__scores_per_site[:, j]
                _weights = np.sum(
                    np.divide(ref_anc_scores_across_sites,
                              total_per_site,
                              out=np.zeros_like(ref_anc_scores_across_sites),
                              where=total_per_site != 0) /
                    Chrom.__samples_per_ref_ancestry_d[self.ref_ancestries[j]])
                normalized_weights[self.ref_ancestries[j]] = _weights

            # add each haplotype's weights
            #interim_res[query_hap_id] = []
            interim_res[query_hap_id] = normalized_weights.copy()

            # reset the scores for next query match evaluation
            Chrom.__scores_per_site.fill(0)
            normalized_weights.clear()

        print("Writing interim files...", flush=True)
        json.dump(interim_res, fw_interm, indent=4)
        fw_interm.close()


def process_per_batch(self, input_folder):
    self.process_all_matches(input_folder)


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
