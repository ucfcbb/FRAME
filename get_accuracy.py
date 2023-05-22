import pandas as pd
import argparse
import os


def calculate_accuracy(inferred_ancestry_file, gt_ancestry_file):
    #gt_file = "../../subsample_indivs_meta.txt"
    df_gt = pd.read_csv(gt_ancestry_file,
                        sep="\t",
                        header=None,
                        names=["Sample", "Ancestry"],
                        index_col=False)
    print(df_gt.head())

    #inferred_file = "./final-inferred-ancestry-allchr-Normed"
    df_inferred = pd.read_csv(inferred_ancestry_file,
                              sep="\t",
                              index_col=False)
    print(df_inferred.head())

    df_merged = pd.merge(df_gt,
                         df_inferred,
                         left_on="Sample",
                         right_on="QuerySample")
    print(df_merged.head())
    print(df_merged[df_merged["Ancestry"] ==
                    df_merged["AssignedAncestry"]].shape[0] / 500)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--i", required=True, help="inferred ancestry file")
    parser.add_argument("--g",
                        required=True,
                        help="ground truth ancestry file")

    args = parser.parse_args()
    inferred_ancestry_file = args.i
    gt_ancestry_file = args.g
    assert os.path.exists(inferred_ancestry_file)
    assert os.path.exists(gt_ancestry_file)
    calculate_accuracy(inferred_ancestry_file, gt_ancestry_file)
