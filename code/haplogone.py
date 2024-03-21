from __future__ import annotations

import os
import pandas as pd
import numpy as np
from cbs.cbs import *
from matplotlib.patches import Rectangle


class VCF:

    def __init__(
        self,
        input_file: str,
        segment_size_threshold: int = 1e6,
        path_to_centromeres_bed: str = "centromeres.bed",
        segmentation_shuffles: int = 1000,
        segmentation_p: float = 0.01,
        validation_shuffles: int = 1000,
        validation_p: float = 0.01,
    ):
        self.input_file = input_file
        self.vcf = None
        self.segment_size_threshold = segment_size_threshold
        self.segmentation_shuffles = segmentation_shuffles
        self.segmentation_p = segmentation_p
        self.validation_shuffles = validation_shuffles
        self.validation_p = validation_p
        self.path_to_centromeres_file = path_to_centromeres_bed
        self.centromeres = self.read_bed(path_to_centromeres_bed)

    def read(self) -> pd.DataFrame:
        """
        Reads a VCF file into a pandas DataFrame.

        Args:
            input_file (str): Path to the input VCF file.

        Returns:
            pd.DataFrame: A DataFrame containing data from the VCF file.
        """

        df_lst = []

        with open(self.input_file) as input_f:
            for current_line in input_f:

                if current_line.startswith("##"):
                    next
                elif current_line.startswith("#CHROM"):
                    colnames = current_line.strip().split("\t")
                else:
                    processed_line_dct = self._process_line(current_line, colnames)
                    df_lst.append(processed_line_dct)

        self.vcf = pd.DataFrame(df_lst)

        return self

    def count_baf(self, ad_cutoff=10) -> pd.DataFrame:
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

        df_vcf = self.vcf

        df_vcf = df_vcf.assign(AD=lambda x: x["AD"].str.split(","))  # split list
        df_vcf["POS"] = pd.to_numeric(df_vcf["POS"])
        df_vcf["DP"] = pd.to_numeric(df_vcf["DP"])
        df_vcf["AD"] = df_vcf.AD.apply(
            lambda x: [num for y in x if (num := int(y)) >= ad_cutoff]
        )  # turn str to int in list

        df_vcf = df_vcf[df_vcf["AD"].map(lambda x: 0 < len(x) < 3)]
        df_vcf["BAF"] = df_vcf.AD.apply(lambda x: max(x) / sum(x))

        self.vcf = df_vcf.reset_index(drop=True)

        return self

    def segment_baf(
        self,
    ) -> pd.DataFrame:
        # TODO write docstring
        df_lst = []

        for chromosome in self.vcf["#CHROM"].unique():

            df_vcf = self.vcf[self.vcf["#CHROM"] == chromosome]

            segments = segment(df_vcf["BAF"].to_numpy())
            segments = validate(df_vcf["BAF"].to_numpy(), segments)  # instability
            segments = self._filter_segments_by_size(segments)

            baf = df_vcf["BAF"].to_numpy()

            if len(segments) > 1:
                baf_segments = self._count_local_mean(segments, baf)

            df_vcf = df_vcf.assign(BAF_segment=baf_segments)

            df_lst.append(df_vcf)

        self.vcf = pd.concat(df_lst, axis=0)

        return self

    def _plot_chromosome(
        self,
        df_vcf: pd.DataFrame,
        name: str,
        chrom_bed,
        chrom_centr: str = None,
        centromeres=None,
    ) -> None:
        """
        TODO write docstring
        """

        fig = plt.figure(figsize=(15, 10))
        fig.suptitle(f"{name}")
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)

        ax1.scatter(
            df_vcf["POS"],
            df_vcf["BAF"],
            s=0.6,
            marker="o",
            label=f"{name}",
        )

        ax1.plot(
            df_vcf["POS"],
            df_vcf["BAF_segment"],
            label=f"{name} segment",
        )

        ax1.plot(
            df_vcf["POS"],
            np.full(len(df_vcf["POS"]), np.mean(df_vcf["BAF_segment"])),
            label=f"{name} mean",
            linestyle="--",
        )
        hist = ax2.hist(
            df_vcf["POS"], bins=180, label=f'"coverage" by POS count"', alpha=0.5
        )
        y, x, _ = hist

        if chrom_centr:
            start_pos = min(
                self.centromeres[self.centromeres["chr"] == chrom_centr]["start"]
            )
            stop_pos = max(
                self.centromeres[self.centromeres["chr"] == chrom_centr]["stop"]
            )

            ax1.add_patch(
                Rectangle(
                    (start_pos, 0.5),
                    stop_pos - start_pos,
                    0.5,
                    facecolor="grey",
                    fill=True,
                    alpha=0.3,
                )
            )

            ax2.add_patch(
                Rectangle(
                    (start_pos, 0),
                    stop_pos - start_pos,
                    max(y),
                    facecolor="grey",
                    fill=True,
                    alpha=0.3,
                )
            )

            ax3.add_patch(
                Rectangle(
                    (start_pos, 0),
                    stop_pos - start_pos,
                    max(df_vcf["DP"]),
                    facecolor="grey",
                    fill=True,
                    alpha=0.3,
                )
            )

        if chrom_bed is not None:
            for start_pos, stop_pos in zip(
                chrom_bed["chromStart"], chrom_bed["chromEnd"]
            ):
                ax1.add_patch(
                    Rectangle(
                        (start_pos, 0.5),
                        stop_pos - start_pos,
                        0.5,
                        facecolor="red",
                        fill=True,
                        alpha=0.3,
                    )
                )

                ax2.add_patch(
                    Rectangle(
                        (start_pos, 0),
                        stop_pos - start_pos,
                        max(y),
                        facecolor="red",
                        fill=True,
                        alpha=0.3,
                    )
                )

                ax3.add_patch(
                    Rectangle(
                        (start_pos, 0),
                        stop_pos - start_pos,
                        max(df_vcf["DP"]),
                        facecolor="red",
                        fill=True,
                        alpha=0.3,
                    )
                )

        ax3.scatter(
            df_vcf["POS"],
            df_vcf["DP"],
            label=f"coverage by DP",
            alpha=0.5,
            s=0.6,
            marker="o",
        )

        ax1.set_xticks(np.arange(0, max(df_vcf["POS"]), step=1e7))
        ax1.set_xticklabels(
            (np.arange(0, max(df_vcf["POS"]), step=1e7) / 1e6).astype(int)
        )
        ax2.set_xticks(np.arange(0, max(df_vcf["POS"]), step=1e7))
        ax2.set_xticklabels(
            (np.arange(0, max(df_vcf["POS"]), step=1e7) / 1e6).astype(int)
        )
        ax3.set_xticks(np.arange(0, max(df_vcf["POS"]), step=1e7))
        ax3.set_xticklabels(
            (np.arange(0, max(df_vcf["POS"]), step=1e7) / 1e6).astype(int)
        )

        fig.legend(loc="lower left")

    def plot_chromosomes(self):
        for chromosome in self.vcf["#CHROM"].unique():
            if chromosome == "chrM":
                continue
            df_vcf = self.vcf[self.vcf["#CHROM"] == chromosome]
            self._plot_chromosome(
                df_vcf,
                name=chromosome,
                chrom_centr=chromosome,
                chrom_bed=self.create_bed(df_vcf),
            )

    def read_bed(self, path_to_bed: str) -> pd.DataFrame:
        """
        Read a BED file and return a pandas DataFrame.

        Args:
            path_to_bed (str): Path to the BED file.

        Returns:
            pd.DataFrame
        """

        df_bed = pd.read_csv(
            path_to_bed,
            sep="\t",
            names=("chr", "start", "stop", "name", "gieStain"),
            skiprows=2,
        )

        return df_bed

    def create_bed(
        self, df_vcf: pd.DataFrame = None, thres: float = 0.9
    ) -> pd.DataFrame:
        """
        Creates a bed pd.DataFrame from a VCF pd.DataFrame with bounds of segments filtered with a given threshold.

        Args:
            df_vcf (pd.DataFrame): Input DataFrame containing VCF data.
            thres (float, optional): Threshold value for BAF_segment. Default is 0.9.

        Returns:
            pd.DataFrame
        """

        if df_vcf is None:
            df_vcf = self.vcf

        df_vcf = df_vcf[df_vcf["BAF_segment"] >= thres]
        df_agg = (
            df_vcf.groupby(["#CHROM", "BAF_segment"])
            .agg({"POS": ["min", "max"], "#CHROM": "first"})["POS"]
            .reset_index()
            .sort_values(["min"])
        ).reset_index(drop=True)

        df_agg["name"] = "LOH"

        df_bed = pd.concat(
            [df_agg["#CHROM"], df_agg["min"], df_agg["max"], df_agg["name"]],
            axis=1,
        )

        df_bed.rename(
            columns={"#CHROM": "#chrom", "min": "chromStart", "max": "chromEnd"},
            inplace=True,
        )

        return df_bed

    def _filter_segments_by_size(
        self,
        segments: np.array,
    ) -> np.array:
        # TODO write docstring

        segment_size_threshold = self.segment_size_threshold

        filtered_segments = [0]
        positions = self.vcf["POS"].to_numpy()

        for indx in range(1, len(segments) - 1):

            if (
                positions[segments[indx]] - positions[filtered_segments[-1]]
                > segment_size_threshold
            ):
                filtered_segments.append(segments[indx])

        filtered_segments.append(segments[-1])

        return filtered_segments

    def _process_line(self, line: str, colnames: list) -> dict:
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

        info_colnames, info_values = [], []

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

    def _count_local_mean(self, segments: np.array, data: np.array) -> np.array:
        """
        Calculates the local mean of BAF (B-allele frequency) values within specified segments.

        Args:
            segments (np.array): An array of segment indices representing the boundaries.
            baf (np.array): An array of BAF values.

        Returns:
            np.array: An array containing the local mean BAF values for each segment.
        """

        local_means_segmented = []

        for indx in range(len(segments) - 1):
            local_mean = data[segments[indx] : segments[indx + 1]].mean()
            local_means_segmented.extend(
                [local_mean] * (segments[indx + 1] - segments[indx])
            )

        return np.array(local_means_segmented)

    def __repr__(self):
        return (
            f"VCF(\n\tinput_file='{self.input_file}',"
            + f"\n\tpath_to_centromeres_file='{self.path_to_centromeres_file}',"
            + f"\n\tsegment_size_threshold={self.segment_size_threshold},"
            + f"\n\tsegmentation_shuffles={self.segmentation_shuffles},"
            + f"\n\tsegmentation_p={self.segmentation_p},"
            + f"\n\tvalidation_shuffles={self.validation_shuffles},"
            + f"\n\tvalidation_p={self.validation_p}\n)"
        )
