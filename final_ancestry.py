import os
import json
#import estimate_final_ancestry
from natsort import natsorted
import collections
import operator

# Number of sites in each chromsome
sitesPerChrom = {}
totalSites = 0
with open('../../numSitesPerChrom', 'r') as f:
    for line in f:
        if 'Chr' in line:
            continue
        toks = line.rstrip().split("\t")
        sitesPerChrom[toks[0]] = int(toks[1])
        totalSites += int(toks[1])

weights_dir = '../weights/'
dirs = os.listdir(weights_dir)
final_ancestry_normed = collections.OrderedDict()

for chrom in natsorted(dirs):
    print(f"Processing {chrom}, num sites = {sitesPerChrom[chrom]}",
          flush=True)
    input_file_path = os.path.join(weights_dir, chrom)
    files = [f for f in os.listdir(input_file_path) if 'faster.json' in f]
    for f in natsorted(files):
        fh = open(os.path.join(input_file_path, f))
        batch_ancestry = json.load(fh)

        for k in batch_ancestry:
            # key_id = indvID without haplotype ID (i.e. -0, -1)
            key_id = k.split('-')[0]

            if key_id not in final_ancestry_normed:
                final_ancestry_normed[key_id] = {}

            ethnicities = batch_ancestry[k]
            for ethnic in ethnicities:
                if ethnic not in final_ancestry_normed[key_id]:
                    final_ancestry_normed[key_id][ethnic] = batch_ancestry[k][
                        ethnic]
                else:
                    final_ancestry_normed[key_id][ethnic] += batch_ancestry[k][
                        ethnic]

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
with open('final-ancestry-normed-2.json', 'w') as fp:
    json.dump(d_descending, fp, indent=4)

# writing to file
#output_file = '/data/pshakya/FindAncestry_withKhazar_pbwtquery/final-inferred-ancestry-allchr-v2'
output_file = 'final-inferred-ancestry-allchr-Normed'
out_write = open(output_file, 'w')
out_write.write('QuerySample\tAssignedAncestry\tCombinedWeight\n')
for topmed_sample in final_ancestry_normed:
    if len(final_ancestry_normed[topmed_sample]) == 0:
        chosen_ethnicity = ("NoMatch", -1)
    else:
        chosen_ethnicity = max(final_ancestry_normed[topmed_sample].items(),
                               key=operator.itemgetter(1))

    out_write.write('{0}\t{1}\t{2:0.4f}\n'.format(topmed_sample,
                                                  chosen_ethnicity[0],
                                                  chosen_ethnicity[1]))
