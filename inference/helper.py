import os
import json
from natsort import natsorted
import collections
import numpy as np


def process_all_weights(num_sites_file, weights_dir):
    sitesPerChrom = {}
    totalSites = 0
    with open(num_sites_file, 'r') as f:
        for line in f:
            if 'Chr' in line:
                continue
            toks = line.rstrip().split("\t")
            sitesPerChrom[toks[0]] = int(toks[1])
            totalSites += int(toks[1])

    weight_files = os.listdir(weights_dir)
    final_ancestry_normed = collections.OrderedDict()

    for weight_file in natsorted(weight_files):
        chrom_num = weight_file.split('-')[0]
        print(f"Processing Chr {chrom_num}, num sites = {sitesPerChrom[chrom_num]}",
            flush=True)

        input_file_path = os.path.join(weights_dir, weight_file)
        fh = open(input_file_path)
        chrom_ancestry = json.load(fh)

        for hap in chrom_ancestry:
            # key_id = indvID without haplotype ID (i.e. -0, -1)
            indv_id = hap.split('-')[0]

            if indv_id not in final_ancestry_normed:
                final_ancestry_normed[indv_id] = {}

            ref_ancestry = chrom_ancestry[hap]
            for anc in ref_ancestry:
                if anc not in final_ancestry_normed[indv_id]:
                    final_ancestry_normed[indv_id][anc] = chrom_ancestry[hap][
                        anc]
                else:
                    final_ancestry_normed[indv_id][anc] += chrom_ancestry[hap][
                        anc]

    # Normalize by the totalNumSites
    for indv_id in final_ancestry_normed:
        for eth in final_ancestry_normed[indv_id]:
            final_ancestry_normed[indv_id][
                eth] = final_ancestry_normed[indv_id][eth] / totalSites

    d_descending = collections.OrderedDict()
    for k in final_ancestry_normed:
        d_descending[k] = dict(
            sorted(final_ancestry_normed[k].items(),
                key=lambda item: item[1],
                reverse=True))
    with open('all-chr-ancestry-normed.json', 'w') as fp:
        json.dump(d_descending, fp, indent=4)


def calculate_proportions(all_chr_ancestry_normed):
    fh = open(all_chr_ancestry_normed, 'r')
    ancestry_scores = json.load(fh)
    ref_populations = []
    for k in ancestry_scores:
        ref_populations = list(ancestry_scores[k].keys())
        break

    print("Total ref populatinos = ", len(ref_populations))
    fw = open('./estimated-ancestry-proportions', 'w')
    fw.write("Sample")  
    for rp in ref_populations:
        fw.write(f"\t{rp}")
    fw.write("\n")

    for indv_id in ancestry_scores:
        fw.write(indv_id)
        total_score = np.sum(list(ancestry_scores[indv_id].values()))
        prop_vals = []
        for hit in ref_populations:
            if total_score != 0:
                proportion_val = ancestry_scores[indv_id][hit] / total_score
            else:
                proportion_val = 0
            prop_vals.append(proportion_val)
            fw.write("\t" + str(proportion_val))

        fw.write("\n")
    fw.close()
