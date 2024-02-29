import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
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


def process_vcf_baf(df_vcf: pd.DataFrame) -> pd.DataFrame:
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
    df_vcf["AD"] = df_vcf.AD.apply(
        lambda x: [num for y in x if (num := int(y)) > 1]
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


def segment_baf(df_vcf: pd.DataFrame, p: float = 0.01) -> pd.DataFrame:
    """
    Segments BAF values in the input DataFrame.

    Args:
        df_vcf (pd.DataFrame): Input DataFrame with a "BAF" column.
        p (float, optional): Threshold for validation. Defaults to 0.01.

    Returns:
        pd.DataFrame: DataFrame with an additional "BAF_segment" column.
    """

    segments = segment(df_vcf["BAF"].to_numpy())
    segments = validate(df_vcf["BAF"].to_numpy(), segments, p=p)  # instability

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

    if any(segment_lengths < threshold):
        pass

    return df_agg
