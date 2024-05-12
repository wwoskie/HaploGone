from __future__ import annotations

import os
import pandas as pd
import random
import string
import numpy as np
import mpld3
import jinja2

from cbs.cbs import *
from datetime import date
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

from pycirclize import Circos


class VCF:
    """
    This is a class for working with the input VCF files. Allows to segment the input data and to vizualize the result.
    """

    def __init__(
        self,
        input_file: str,
        segment_size_threshold: int = 1e6,
        path_to_centromeres_bed: str = "centromeres.bed",
        segmentation_shuffles: int = 1000,
        segmentation_p: float = 0.01,
        validation_shuffles: int = 1000,
        validation_p: float = 0.01,
        baf_freq_threshold: float = 0.95,
        output_dir=None,
    ):
        self.input_file = input_file
        self.vcf = None
        self.bed = None
        self.segment_size_threshold = segment_size_threshold
        self.segmentation_shuffles = segmentation_shuffles
        self.segmentation_p = segmentation_p
        self.validation_shuffles = validation_shuffles
        self.validation_p = validation_p
        self.baf_freq_threshold = baf_freq_threshold
        self.path_to_centromeres_file = path_to_centromeres_bed
        self.centromeres = self.read_bed(path_to_centromeres_bed)

        if output_dir is None:
            current_date = date.today().strftime("%Y-%m-%d")
            base = os.path.basename(input_file)
            file = os.path.splitext(base)[0]
            rand_id = ''.join(random.choice(
                string.ascii_uppercase + string.digits) for _ in range(6))

            output_dir = f"{current_date}_{file}_{rand_id}"

        self.output_dir = output_dir

    def read(self) -> pd.DataFrame:
        """
        Converts a VCF file into a pandas DataFrame.

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
                    processed_line_dct = self._process_line(
                        current_line, colnames)
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

        df_vcf = df_vcf.assign(
            AD=lambda x: x["AD"].str.split(","))  # split list
        df_vcf["POS"] = pd.to_numeric(df_vcf["POS"])
        df_vcf["DP"] = pd.to_numeric(df_vcf["DP"])
        df_vcf["AD"] = df_vcf.AD.apply(
            lambda x: [num for y in x if (num := int(y)) >= ad_cutoff]
        )  # turn str to int in list

        df_vcf = df_vcf[df_vcf["AD"].map(lambda x: 0 < len(x) < 3)]
        df_vcf["BAF"] = df_vcf.AD.apply(lambda x: max(x) / sum(x))

        self.vcf = df_vcf.reset_index(drop=True)

        return self

    def segment_baf(self) -> pd.DataFrame:
        """
        Segments BAF values in the input DataFrame.

        Args:
            df_vcf (pd.DataFrame): Input DataFrame with a "BAF" column.

        Returns:
            pd.DataFrame: The DataFrame with an additional "BAF_segment" column.
        """
        df_lst = []

        for chromosome in self.vcf["#CHROM"].unique():

            df_vcf = self.vcf[self.vcf["#CHROM"] == chromosome]
            baf = df_vcf["BAF"].to_numpy()

            segments = segment(baf)
            segments = validate(baf, segments)  # instability
            segments = self._filter_segments_by_size(segments)

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
        save_plot = False,
        to_html=False,
    ) -> None:
        """
        Vizuales the result of the segmentation.

        Args: 
            df_vcf (pd.DataFrame): Input DataFrame with a "BAF_segment" column.
        
        Return:
            The plot of the chromosome data with segmentation.
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

        if chrom_centr:
            start_pos = min(
                self.centromeres[self.centromeres["chr"]
                                 == chrom_centr]["start"]
            )
            stop_pos = max(
                self.centromeres[self.centromeres["chr"]
                                 == chrom_centr]["stop"]
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

        ax3.scatter(
            df_vcf["POS"],
            df_vcf["DP"],
            label=f"coverage by DP",
            alpha=0.5,
            s=0.8,
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


        if save_plot:
            self._check_and_create_output_dir()
            fig.savefig(f'{self.output_dir}/{name}.png')   # save the figure to file
            plt.close(fig)
            return f'{name}.png'
        
        if to_html:
            plt.close(fig) 
            return mpld3.fig_to_html(fig)

    def plot_chromosomes(self, save_plot=False, to_html=False,):
        """
        Vizualizes the whole result of the segmentation for each chromosome in the dataset.
        
        Args: 
            df_vcf (pd.DataFrame): Input DataFrame with a "BAF_segment" column.
        
        Return:
            The plots of the chromosome data with segmentation.

        """
        if self.bed is None:
            self.bed = self.create_bed(self.vcf)

        html_list = []
        plot_paths = []


        for chromosome in self.vcf["#CHROM"].unique():
            if chromosome == "chrM":
                continue
            df_vcf = self.vcf[self.vcf["#CHROM"] == chromosome]
            

            if to_html:
                html_list.append(
                    self._plot_chromosome(
                    df_vcf,
                    name=chromosome,
                    chrom_centr=chromosome,
                    chrom_bed=self.bed[self.bed["#chrom"] == chromosome],
                    save_plot=save_plot,
                    to_html=to_html,
                    )
                )

            elif save_plot:
                plot_paths.append(
                    self._plot_chromosome(
                    df_vcf,
                    name=chromosome,
                    chrom_centr=chromosome,
                    chrom_bed=self.bed[self.bed["#chrom"] == chromosome],
                    save_plot=save_plot,
                )
            )

        if to_html:
            return html_list

        if save_plot:
            return plot_paths


    def plot_chromosomes_circular(self, save_plot=False, to_html=False):
        sectors = {}

        for chrom in self.vcf["#CHROM"].unique():
            sectors[chrom] = max(self.vcf[self.vcf["#CHROM"] == chrom]["POS"])

        if "chrM" in sectors:
            del sectors["chrM"]

        from pycirclize import Circos

        circos = Circos(
            sectors,
            space=1,
            endspace=False,
        )

        for sector in circos.sectors:
            chrom_bed = self.bed[self.bed["#chrom"] == sector.name]
            sector.axis(lw=0)
            sector.text(f"{sector.name.replace("chr", "")}", size=8)

            track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
            track1.axis()

            track1.scatter((self.vcf[self.vcf["#CHROM"] == sector.name]["POS"]).to_numpy(
            ), (self.vcf[self.vcf["#CHROM"] == sector.name]["BAF"]).to_numpy(), s=0.1, color="violet", alpha=0.3)
            track1.line((self.vcf[self.vcf["#CHROM"] == sector.name]["POS"]).to_numpy(
            ), (self.vcf[self.vcf["#CHROM"] == sector.name]["BAF_segment"]).to_numpy(), color="blue")

            for start_pos, stop_pos in zip(chrom_bed["chromStart"], chrom_bed["chromEnd"]):
                sector.rect(start=start_pos, end=stop_pos, r_lim=(80, 100), color="red", alpha=0.3)

            # add centromeres
            start_pos = min(
                self.centromeres[self.centromeres["chr"]
                                 == sector.name]["start"]
            )
            stop_pos = max(
                self.centromeres[self.centromeres["chr"]
                                 == sector.name]["stop"]
            )

            sector.rect(start=start_pos, end=stop_pos, r_lim=(80, 100), color="grey", alpha=0.7)

        fig = circos.plotfig()

        legend_handles = [
            Line2D([], [], color="violet", marker="o", label="BAF", ms=3, ls="None"),
            Line2D([], [], color="blue", label="BAF segments"),
            Patch(color="red", label="LOH"), 
            Patch(color="grey", label="centromere")
        ]

        _ = circos.ax.legend(
            handles=legend_handles,
            bbox_to_anchor=(0.5, 0.5),
            loc="center",
            fontsize=8,
        )

        

        if save_plot:
            self._check_and_create_output_dir()
            fig.savefig(f'{self.output_dir}/circular_plot.png', dpi=300)   # save the figure to file
            plt.close(fig)

        if to_html:
            return mpld3.fig_to_html(fig)

    def generate_report(self):
        self.plot_chromosomes_circular()
        chrom_names = self.vcf["#CHROM"].unique()
        circle_plot_path = f'circular_plot.png'
        plots = self.plot_chromosomes(save_plot=True)

        templateLoader = jinja2.FileSystemLoader(searchpath="./")
        templateEnv = jinja2.Environment(loader=templateLoader)
        TEMPLATE_FILE = "index.html"
        template = templateEnv.get_template(TEMPLATE_FILE)
        output = template.render(circle_plot_path = circle_plot_path, 
        plots = plots, chrom_names = chrom_names)

        with open(f'{self.output_dir}/report.html', 'w') as f:
            f.write(output)

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
        self, df_vcf: pd.DataFrame = None, thres: float = None
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

        if thres is None:
            thres = self.baf_freq_threshold

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
            columns={"#CHROM": "#chrom",
                     "min": "chromStart", "max": "chromEnd"},
            inplace=True,
        )
        df_bed = df_bed.sort_values(by=["#chrom"]).reset_index(drop=True)

        if self.bed is None:
            self.bed = df_bed

        return self

    def save_bed(self):
        if self.bed is None:
            self.create_bed()
        self._check_and_create_output_dir()
        self.bed.to_csv(f'{self.output_dir}/output.bed', index=False, sep='\t')  

        return self


    def _filter_segments_by_size(self, segments: np.array) -> np.array:
        """
        Filters the segments by size. The segments smaller than the threshold would be removed.

        Args:
            A column with number of  variant position.
        Return:
            A list of filtered segments.
        """

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
            local_mean = data[segments[indx]: segments[indx + 1]].mean()
            local_means_segmented.extend(
                [local_mean] * (segments[indx + 1] - segments[indx])
            )

        return np.array(local_means_segmented)

    def _check_and_create_output_dir(self):
        if not os.path.isdir(self.output_dir):
                os.mkdir(self.output_dir)


    def __repr__(self):
        return (
            f"VCF(\n\tinput_file='{self.input_file}',"
            + f"\n\tpath_to_centromeres_file='{self.path_to_centromeres_file}',"
            + f"\n\tbaf_freq_threshold={self.baf_freq_threshold},"
            + f"\n\tsegment_size_threshold={self.segment_size_threshold},"
            + f"\n\tsegmentation_shuffles={self.segmentation_shuffles},"
            + f"\n\tsegmentation_p={self.segmentation_p},"
            + f"\n\tvalidation_shuffles={self.validation_shuffles},"
            + f"\n\tvalidation_p={self.validation_p}\n)"
        )
