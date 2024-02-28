import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


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


def process_vcf_baf(table):
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

    table = table.assign(AD=lambda x: x["AD"].str.split(","))  # split list
    table["POS"] = pd.to_numeric(table["POS"])
    table["AD"] = table.AD.apply(
        lambda x: [int(y) for y in x]
    )  # turn str to int in list
    table["AD"] = table.AD.apply(
        lambda x: [y for y in x if y != 0 and y != 1]
    )  # filter zeros (and ones) from list
    table = table[
        (table["AD"].map(lambda x: len(x)) > 0)
        & (table["AD"].map(lambda x: len(x)) < 3)
    ]
    table["BAF"] = table.AD.apply(lambda x: x[0] / sum(x))

    return table
