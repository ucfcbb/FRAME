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

    def __init__(self, ref_popln_file, query_sample_file, output_folder,  numSiteFile, chrom_num):
        self.chrom_num = chrom_num
        print(f"---------- Initializing CHR {self.chrom_num} -----------",
              flush=True)

        # Load haplotype to Individual mapping for query & ref. samples
        self.query_samples = self.get_query_samples(query_sample_file)

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

    def get_query_samples(self, query_sample_file):
        '''
            returns the mapping of query haplotype index to query individual id
        '''
        query_samples = []
        with open(query_sample_file, 'r') as f:
            for line in f:
                toks = line.rstrip().split()
                idd = toks[0]
                query_samples.append(idd+'-0')
                query_samples.append(idd+'-1')
        return query_samples

    def load_all_required_meta_info(self, ref_popln_file, delim="\t"):
        print("Loading reference popuation...", flush=True)
        idx = 0
        with open(ref_popln_file, 'r') as f:
            for line in f:
                if line.startswith('Sample'):
                    continue
                toks = line.rstrip().split(delim)
                sample_id = toks[0]
                population = toks[1]
                Chrom.ref_ancestry_of_sample_d[sample_id] = population
                Chrom.samples_per_ref_ancestry_d[population] += 1

                # assign some arbritrary index for each ref-popln
                # to initilize and index the scores_per_site array
                if population not in Chrom.idx_of_ref_ancestry_d:
                    Chrom.idx_of_ref_ancestry_d[population] = idx
                    idx += 1

    def process_all_matches(self, chrom_IBD_file):
        ''' 
            Process all the matching segments for a chromosome
        '''
        # create the interm output file
        fw_interm = open(self.interm_output_path, 'w')

        # initialize interim_res
        interim_res = collections.OrderedDict()
        for qs in self.query_samples:
            interim_res[qs] = collections.OrderedDict()
            for ra in self.ref_ancestries:
                interim_res[qs][ra] = 0.0

        # Update #hits per query haplotype
        prev_query_sample = ''
        with open(chrom_IBD_file, 'r') as f:
            for line in f:
                toks = line.rstrip().split()
                query_sample = toks[0]
                ref_sample, ref_hap_id = toks[1].split("-")
                start_site = int(toks[2])
                end_site = int(toks[3])  # exclusive

                if prev_query_sample == '':
                    prev_query_sample = query_sample
                elif prev_query_sample != query_sample:
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

                        # add each haplotype's weights
                        interim_res[prev_query_sample][self.ref_ancestries[j]] = weights
                    
                    # reset the scores for next query match evaluation
                    Chrom.scores_matrix.fill(0)
                    prev_query_sample = query_sample
                    
                # update for the current query haplotype
                array_idx_ref_ancestry = Chrom.idx_of_ref_ancestry_d[Chrom.ref_ancestry_of_sample_d[ref_sample]]
                Chrom.scores_matrix[start_site:end_site,
                                        array_idx_ref_ancestry] += 1
            
        # handle for the boundary case
        # last query sample(s) is(are) the same
        total_per_site = np.sum(Chrom.scores_matrix, axis=1)
        for j in range(self.total_num_uniq_ref_ancestries):
            ref_ancestry_scores_across_sites = Chrom.scores_matrix[:, j]
            weights = np.sum(
                np.divide(ref_ancestry_scores_across_sites,
                            total_per_site,
                            out=np.zeros_like(ref_ancestry_scores_across_sites),
                            where=total_per_site != 0) /
                Chrom.samples_per_ref_ancestry_d[self.ref_ancestries[j]])
            interim_res[query_sample][self.ref_ancestries[j]] = weights

        print("Writing interim files...", flush=True)
        json.dump(interim_res, fw_interm, indent=4)
        fw_interm.close()
