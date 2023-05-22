import os
import json
import numpy as np

final_ancestry_normed = "./final-ancestry-normed-2.json"
fh = open(final_ancestry_normed, 'r')
ancestry_scores = json.load(fh)

fw = open('./estimated-ancestry-proportions', 'w')
fw.write("Sample\tEUR\tAFR\tAMR\tEAS\tSAS\n")

ref_ancestires_order = ["EUR", "AFR", "AMR", "EAS", "SAS"]
for indv_id in ancestry_scores:
    fw.write(indv_id )
    total_score = np.sum(list(ancestry_scores[indv_id].values()))
    for hit in ref_ancestires_order:
        proportion_val = ancestry_scores[indv_id][hit] / total_score
        fw.write("\t" + str(proportion_val))
    fw.write("\n")
fw.close()
