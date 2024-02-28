import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cbs


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


def filter_vcf(input_path: str) -> pd.DataFrame:
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
