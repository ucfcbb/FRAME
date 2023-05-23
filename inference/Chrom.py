import os
import collections
import numpy as np
from natsort import natsorted
import argparse
import json


class Chrom:
    # store all the calculated ancestral scores
    # 2D numpy array : #sitesInChrom X #UniqRefAncestries
    scores_matrix = None

    # Reference metadata information #
    ref_ancestry_of_sample_d = collections.defaultdict()
    idx_of_ref_ancestry_d = collections.OrderedDict()
    samples_per_ref_ancestry_d = collections.defaultdict(lambda: 0)
    num_sites_per_chrom = collections.defaultdict()

    def __init__(self, chrom_IBD_file, output_folder, ref_popln_file, numSiteFile, chrom_num):
        self.chrom_num = chrom_num
        print(f"------------- Initializing {self.chrom_num} -----------",
              flush=True)

        # Load haplotype to Individual mapping for query & ref. samples
        #self.haps2indivs_query = self.get_query_indivs()
        self.haps2indivs_ref, self.uniq_indivs = self.get_indivs_ref(ref_popln_file)

        # loading meta information
        self.get_num_sites_per_chrom(numSiteFile)
        self.load_all_required_meta_info(ref_popln_file)

        self.total_num_sites = Chrom.num_sites_per_chrom[self.chrom_num]
        self.total_num_uniq_ref_ancestries = len(Chrom.idx_of_ref_ancestry_d)
        self.ref_ancestries = list(Chrom.idx_of_ref_ancestry_d.keys())

        # create 2D numpy array : #sitesInChrom X #UniqRefAncestries
        Chrom.scores_matrix = np.zeros(
            shape=(self.total_num_sites, self.total_num_uniq_ref_ancestries))

        # setup output directory and output weights file
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        self.interm_output_path = os.path.join(
            output_folder,
            self.chrom_num + '-interm-ancestry-result-faster.json')


    def get_num_sites_per_chrom(self, numSiteFile):
        ''' 
            Stores the # of sites present in every chromosome's VCF file in a dict.
        '''
        with open(numSiteFile, 'r') as f:
            for line in f:
                if 'Chr' in line:
                    continue
                toks = line.strip().split()
                Chrom.num_sites_per_chrom[toks[0]] = int(toks[1])

    def get_query_hap_idx(self, input_file):
        '''
            returns the query haplotype's index
        '''
        return input_file.split('_')[-1]

    def get_query_indivs(self, ):
        '''
            returns the mapping of query haplotype index to query individual id
        '''
        hap_indices = {}
        hap_id = 0
        prefix = 'indivs-list'
        with open(prefix, 'r') as f:
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
        print("Loading reference popuation...", flush=True)
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

                Chrom.ref_ancestry_of_sample_d[sample_id] = population

                # store # of samples per ref_popln
                Chrom.samples_per_ref_ancestry_d[population] += 1

                # assign some arbritrary index for each ref-popln
                # to initilize and index the scores_per_site array
                if population not in Chrom.idx_of_ref_ancestry_d:
                    Chrom.idx_of_ref_ancestry_d[population] = idx
                    idx += 1

    def process_all_matches(self, chrom_IBD_file):
        ''' 
            Process all the matching segments for a chromosome (w/o batching)
        '''

        # create the interm output file
        fw_interm = open(self.interm_output_path, 'w')

        # iterate through all match files
        interim_res = {}
        normalized_weights = {}
        #query_hap_idx = self.get_query_hap_idx(chrom_IBD_file)
        #query_hap_id = self.haps2indivs_query[query_hap_idx]

        # Update #hits per site
        min_site = -1
        max_site = -1
        with open(chrom_IBD_file, 'r') as f:
            for line in f:
                toks = line.rstrip().split()
                query_sample, query_hap_id = toks[0].split("-")
                ref_sample, ref_hap_id = toks[1].split("-")
                start_site = int(toks[2])
                end_site = int(toks[3])  # exclusive

                if min_site < 0:
                    if start_site > min_site:
                        min_site = start_site
                else:
                    if start_site < min_site:
                        min_site = start_site

                if end_site > max_site:
                    max_site = end_site

                #hit_ref_hap_idx = self.haps2indivs_ref[ref_hap_idx]
                #hit_ref_hap_sample = hit_ref_hap_idx[:-2]

                array_idx_ref_ancestry = Chrom.idx_of_ref_ancestry_d[ref_sample]
                Chrom.scores_matrix[start_site:end_site,
                                        array_idx_ref_ancestry] += 1

        # process weights
        total_per_site = np.sum(Chrom.scores_matrix, axis=1)
        for j in range(self.total_num_uniq_ref_ancestries):
            ref_ancestry_scores_across_sites = Chrom.scores_matrix[:, j]
            weights = np.sum(
                np.divide(ref_ancestry_scores_across_sites,
                            total_per_site,
                            out=np.zeros_like(ref_ancestry_scores_across_sites),
                            where=total_per_site != 0) /
                Chrom.samples_per_ref_ancestry_d[self.ref_ancestries[j]])

            normalized_weights[self.ref_ancestries[j]] = weights

        # add each haplotype's weights
        interim_res[query_hap_id] = normalized_weights.copy()

        # reset the scores for next query match evaluation
        Chrom.scores_matrix.fill(0)
        normalized_weights.clear()

        print("Writing interim files...", flush=True)
        json.dump(interim_res, fw_interm, indent=4)
        fw_interm.close()