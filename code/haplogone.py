import os
import pandas as pd
import numpy as np
from cbs.cbs import *


def process_line(line: str, colnames: list) -> dict:
    """
    Processes a single line from a VCF file.

    Args:
        line (str): A single line from the VCF file.
        colnames (list): List of column names extracted from the header.

    Returns:
        dict: A dictionary containing processed information from the line.
            The keys correspond to column names, and the values are the
            corresponding values from the line.
    """

    line_values = line.strip().split("\t")

    sample_values = line_values.pop().split(":")
    format_colnames = line_values.pop().split(":")
    info = line_values.pop().split(";")

    info_colnames = []
    info_values = []

    for element in info:
        if "=" in element:
            info_colname, info_value = element.split("=")
        else:
            info_colname, info_value = element, ""

        info_colnames.append(info_colname)
        info_values.append(info_value)

    colnames = colnames[:-3]
    colnames.extend(info_colnames + format_colnames)
    line_values.extend(info_values + sample_values)
    output_dct = dict(zip(colnames, line_values))

    return output_dct


def read_vcf(input_path: str) -> pd.DataFrame:
    """
    Reads a VCF file and filters relevant data into a pandas DataFrame.

    Args:
        input_path (str): Path to the input VCF file.

    Returns:
        pd.DataFrame: A DataFrame containing processed data from the VCF file.
    """

    df_lst = []

    with open(input_path) as input_f:
        for current_line in input_f:

            if current_line.startswith("##"):
                next
            elif current_line.startswith("#CHROM"):
                colnames = current_line.strip().split("\t")
            else:
                processed_line_dct = process_line(current_line, colnames)
                df_lst.append(processed_line_dct)

    return pd.DataFrame(df_lst)


def process_vcf_baf(df_vcf: pd.DataFrame, ad_cutoff=0) -> pd.DataFrame:
    """
    Processes variant data from a VCF (Variant Call Format) DataFrame to calculate B-Allele Frequency (BAF).

    Args:
        df_vcf (pd.DataFrame): A DataFrame containing VCF variant data with the following columns:
            - 'AD': Allelic depths (comma-separated string of integers).
            - 'POS': Variant position (numeric).

    Returns:
        pd.DataFrame: Processed DataFrame with the following additional column:
            - 'BAF': B-Allele Frequency calculated as the maximum allelic depth divided by the total allelic depth.
    """

    df_vcf = df_vcf.assign(AD=lambda x: x["AD"].str.split(","))  # split list
    df_vcf["POS"] = pd.to_numeric(df_vcf["POS"])
    df_vcf["DP"] = pd.to_numeric(df_vcf["DP"])
    df_vcf["AD"] = df_vcf.AD.apply(
        lambda x: [num for y in x if (num := int(y)) >= ad_cutoff]
    )  # turn str to int in list

    df_vcf = df_vcf[df_vcf["AD"].map(lambda x: 0 < len(x) < 3)]
    df_vcf["BAF"] = df_vcf.AD.apply(lambda x: max(x) / sum(x))

    return df_vcf.reset_index(drop=True)


def generate_baf_segments(segments: np.array, baf: np.array) -> np.array:
    """
    Calculates the local mean of BAF (B-allele frequency) values within specified segments.

    Args:
        segments (np.array): An array of segment indices representing the boundaries.
        baf (np.array): An array of BAF values.

    Returns:
        np.array: An array containing the local mean BAF values for each segment.
    """

    baf_segments = []

    for indx in range(len(segments) - 1):
        local_mean = baf[segments[indx] : segments[indx + 1]].mean()
        baf_segments.extend([local_mean] * (segments[indx + 1] - segments[indx]))

    return np.array(baf_segments)


def segment_baf(
    df_vcf: pd.DataFrame, shuffles: int = 1000, p: float = 0.01
) -> pd.DataFrame:
    """
    Segments BAF values in the input DataFrame.

    Args:
        df_vcf (pd.DataFrame): Input DataFrame with a "BAF" column.
        p (float, optional): Threshold for validation. Defaults to 0.01.

    Returns:
        pd.DataFrame: DataFrame with an additional "BAF_segment" column.
    """

    segments = segment(df_vcf["BAF"].to_numpy(), shuffles=shuffles, p=p)
    segments = validate(
        df_vcf["BAF"].to_numpy(), segments, shuffles=shuffles, p=p
    )  # instability

    baf = df_vcf["BAF"].to_numpy()

    if len(segments) > 1:
        baf_segments = generate_baf_segments(segments, baf)

    df_vcf = df_vcf.assign(BAF_segment=baf_segments)

    return df_vcf


def filter_segments_by_size(df_vcf: pd.DataFrame, threshold: int) -> pd.DataFrame:

    df_agg = (
        df_vcf.groupby("BAF_segment")
        .agg({"POS": ["min", "max"]})["POS"]
        .reset_index()
        .sort_values(["min"])
    )

    segments = np.array(df_agg["BAF_segment"])
    segment_lengths = np.array(df_agg["max"] - df_agg["min"])

    df_agg["segment_length"] = segment_lengths

    new_maxes = []
    new_mins = []
    new_length = []

    if any(segment_lengths < threshold):
        for indx in range(len(df_agg["segment_length"])):
            min_ = df_agg["min"][indx]
            max_ = df_agg["max"][indx]
            length = df_agg["segment_length"][indx]

            if length > threshold:
                new_mins.append()
                new_maxes.append

    return df_agg


def read_and_segment(
    path: str,
    ad_cutoff: int = 10,
    shuffles: int = 1000,
    p: float = 1e-5,
    chrom_subset: str = None,
) -> pd.DataFrame:
    """
    Reads a VCF file from the specified path, processes it, and segments the data.

    Args:
        path (str): Path to the VCF file.
        ad_cutoff (int, optional): Minimum allele depth cutoff. Defaults to 10.
        shuffles (int, optional): Number of shuffles for BAF segmentation. Defaults to 1000.
        p (float, optional): Significance threshold for segmentation. Defaults to 1e-5.
        chrom_subset (str, optional): Subset of chromosomes to consider. Defaults to None.

    Returns:
        pd.DataFrame: Processed and segmented VCF data.
    """

    vcf_df = read_vcf(path)
    if chrom_subset:
        vcf_df = vcf_df[vcf_df["#CHROM"] == chrom_subset]
    vcf_df = process_vcf_baf(vcf_df, ad_cutoff)
    vcf_df = segment_baf(vcf_df, shuffles, p)

    return vcf_df


def plot_chromosome(df_vcf: pd.DataFrame, name: str) -> None:
    """
    Plots B-Allele Frequency (BAF) data for a given chromosome.

    Args:
        df_vcf (pd.DataFrame): DataFrame containing VCF data.
        name (str): Name of the chromosome.

    Returns:
        None: Displays the BAF plot.
    """

    fig = plt.figure(figsize=(20, 15))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax1.scatter(
        df_vcf["POS"],
        df_vcf["BAF"],
        s=0.6,
        # c="b",
        marker="o",
        label=f"{name}",
    )

    ax1.plot(
        df_vcf["POS"],
        df_vcf["BAF_segment"],
        # c="black",
        label=f"{name} segment",
    )

    ax1.plot(
        df_vcf["POS"],
        np.full(len(df_vcf["POS"]), np.mean(df_vcf["BAF_segment"])),
        label=f"{name} mean",
        linestyle="--",
    )

    ax2.hist(df_vcf["POS"], bins=180, label=f'"coverage" by POS count"', alpha=0.5)
    ax3.scatter(df_vcf["POS"], df_vcf["DP"], label=f"coverage by DP", alpha=0.5)

    ax1.set_xticks(np.arange(0, max(df_vcf["POS"]), step=1e7))
    ax2.set_xticks(np.arange(0, max(df_vcf["POS"]), step=1e7))
    ax3.set_xticks(np.arange(0, max(df_vcf["POS"]), step=1e7))

    fig.legend(loc="lower left")
