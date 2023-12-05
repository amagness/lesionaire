import sys

import argparse
import glob
import os

import alphashape
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from descartes import PolygonPatch
from sklearn.cluster import DBSCAN
from tqdm import *

sys.path.append("./")
from lesionaire.plots import plot_clusters_pub_proportional


class lesionData:
    """
    The lesionaire lesionData class.
    """

    def __init__(
        self,
        data,
        x_id_col,
        y_id_col,
        clustering_id_col,
        class_id,
        cluster_eps=35,
        min_s=1,
        cluster_alpha=0.05,
        tissue_alpha=0.05,
        image_id_col='Image',
        lobe_id_col=None,
        cell_area_col=None,
        lesion_threshold=5,
        area_scale=1e6,
    ) -> None:
        """Initialise lesionData object."""

        self.data = data
        self.image_id_col = image_id_col
        self.lobe_id_col = lobe_id_col
        self.x_id_col = x_id_col
        self.y_id_col = y_id_col
        self.clustering_id_col = clustering_id_col
        self.class_id = class_id
        self.cluster_eps = cluster_eps
        self.cluster_alpha = cluster_alpha
        self.tissue_alpha = tissue_alpha
        self.min_s = min_s  # minimum number of cells in a cluster for DBSCAN clustering
        self.lesion_threshold = lesion_threshold  #
        self.lesion_plot_cutoff = (
            self.lesion_threshold
        )  # clusters with less than this number of cells will be plotted as gray points
        self.cell_area_col = cell_area_col
        self.area_scale = area_scale  # scale factor to convert area units to desired units e.g. 1e6 to convert from lesions/µm^2 to lesions/mm^2

        self.find_lesions()
        print("Found {} lesions.".format(self.total_lesions))

    def find_lesions(self):
        """Find lesions with density-based clustering (DBSCAN).

        Returns:
            pd.DataFrame: Original data with lesion ID column added.
        """
        self.lesions = process(
            self.data,
            self.image_id_col,
            self.x_id_col,
            self.y_id_col,
            self.clustering_id_col,
            cluster_class=self.class_id,
            cluster_eps=self.cluster_eps,
            min_s=self.min_s,
        )
        self.total_lesions = self.lesions["lesion_id"].nunique()
        self.total_clustered_cells = len(self.lesions)

    def get_lesions(self):
        return self.lesions

    def cell_counts_per_lesion(self):
        if self.lobe_id_col:
            self.counts = (
                self.lesions.groupby([self.image_id_col, self.lobe_id_col, "lesion_id"])
                .size()
                .reset_index(name="cell_count")
            )
        else:
            self.counts = (
                self.lesions.groupby([self.image_id_col, "lesion_id"])
                .size()
                .reset_index(name="cell_count")
            )
        return self.counts

    # def count_lesions(self):

    def summary(self):
        """
        Returns a summary of the lesions in the dataset.
        """

        self.summary = (
            self.cell_counts_per_lesion()
            .groupby(self.image_id_col)
            .agg({"cell_count": ["mean", "std"]})
        )
        self.summary.columns = [
            "_".join(col).strip() for col in self.summary.columns.values
        ]
        self.summary = self.summary.reset_index()
        rename_dict = dict(
            {
                "cell_count_mean": "mean_cells_per_lesion",
                "cell_count_std": "std_cells_per_lesion",
            }
        )
        self.summary.rename(columns=rename_dict, inplace=True)

        # check if class has tissue area attribute:
        if not hasattr(self, "tissue_area"):
            self.get_tissue_boundary()

        self.summary["total_lesions"] = self.lesions["lesion_id"].nunique()
        self.summary["tissue_area"] = self.tissue_area
        self.summary["lesions_per_unit_tissue_area"] = (
            self.summary["total_lesions"][0] / self.summary["tissue_area"][0]
        ) * self.area_scale
        if self.cell_area_col:  # calculate cell area from single cell area measurements if cell area column is provided
            self.summary["total_cell_area"] = self.data[self.cell_area_col].sum()
            self.summary["lesions_per_unit_cell_area"] = (
                self.summary["total_lesions"][0] / self.summary["total_cell_area"][0]
            ) * self.area_scale
            self.total_cell_area = self.summary["total_cell_area"][0]
        self.summary["total_clustered_cells"] = self.total_clustered_cells
        self.summary["min_lesion_size"] = self.min_s
        self.summary["max_lesion_size"] = self.lesions["Cluster N cells"].max()
        self.summary["cluster_eps"] = self.cluster_eps

        return self.summary

    def get_tissue_boundary(self):
        """
        Function to calculate the area of a set of points with the alpha shape algorithm.
        """

        # create alphashape:
        print("Calculating alphashape tissue boundary...")
        points = self.data[[self.x_id_col, self.y_id_col]].values
        self.tissue_boundary = alphashape.alphashape(points, self.tissue_alpha)
        self.tissue_area = self.tissue_boundary.area
        print("Done.")
        return self.tissue_boundary

    def plot_lesions(self, outdir=None, sample_name=""):
        """
        Plots identified lesions with an alphashape pertaining to each lesion over the background tissue.

        """

        sns.set_style("white")

        ## make alphashape and plot outline of lung
        if not hasattr(self, "tissue_boundary"):
            self.get_tissue_boundary()

        image_shape = (self.data[self.y_id_col].max(), self.data[self.x_id_col].max())

        # create figure:
        scale = 5e-3
        fig, ax = plt.subplots(
            nrows=1, ncols=1, figsize=(image_shape[1] * scale, image_shape[0] * scale)
        )

        # plot lung:
        ax.add_patch(PolygonPatch(self.tissue_boundary, alpha=0.3)).set(
            facecolor="#f0b9e7", edgecolor="k"
        )

        ## plot small clusters (<5 cutoff) as grey:
        small_clusters = self.lesions[
            self.lesions["Cluster N cells"] < self.lesion_plot_cutoff
        ]
        ax.scatter(
            small_clusters[self.x_id_col].values,
            small_clusters[self.y_id_col].values,
            alpha=0.1,
            color="k",
            s=small_clusters["Cluster N cells"].values * 250,
        )
        ax.scatter(
            small_clusters[self.x_id_col].values,
            small_clusters[self.y_id_col].values,
            alpha=1,
            s=1,
            color="k",
        )

        large_clusters = self.lesions[
            self.lesions["Cluster N cells"] >= self.lesion_plot_cutoff
        ]

        # get unique cluster labels:
        unique_cluster_labels = large_clusters["Cluster ID"].unique()

        # loop through labels:
        for label in unique_cluster_labels:
            # get df for each unique label:
            cluster_df = large_clusters[large_clusters["lesion_id"] == label]

            # get points of cluster:
            points = cluster_df[[self.x_id_col, self.y_id_col]].values

            # only proceed if cells exist:
            if len(points.tolist()) > 0:
                if label == -1:
                    ax.scatter(
                        cluster_df[self.x_id_col].values,
                        cluster_df[self.y_id_col].values,
                        alpha=0.75,
                        color="k",
                        s=50,
                    )

                else:
                    # create alphashape:
                    alpha_shape = alphashape.alphashape(points, self.cluster_alpha)

                    if (
                        alpha_shape.geom_type in ["Polygon", "MultiPolygon"]
                    ):  # only process Polygon and Multipolygons i.e. ignore lines of cells which cannot contain other cells
                        # plot points and add patch of alphashape:
                        ax.scatter(
                            cluster_df[self.x_id_col].values,
                            cluster_df[self.y_id_col].values,
                            s=cluster_df["Cluster N cells"].values * 250,
                            alpha=0.05,
                        )
                        ax.scatter(
                            cluster_df[self.x_id_col].values,
                            cluster_df[self.y_id_col].values,
                            alpha=1,
                            s=1,
                            color="k",
                        )
                        ax.add_patch(PolygonPatch(alpha_shape, alpha=0.3)).set(
                            facecolor="k", edgecolor="w"
                        )
                    else:
                        ax.scatter(
                            cluster_df[self.x_id_col].values,
                            cluster_df[self.y_id_col].values,
                            alpha=0.1,
                            color="k",
                            s=150,
                        )
                        ax.scatter(
                            cluster_df[self.x_id_col].values,
                            cluster_df[self.y_id_col].values,
                            alpha=1,
                            color="k",
                            s=1,
                        )

        # update plot with title etc and save:
        # title = '{}_{}_alpha_{}.png'.format(self.sample_name, self.class_id, self.cluster_alpha)
        ax.set_xlabel("Cell position X µm", fontsize=75)
        ax.set_ylabel("Cell position Y µm", fontsize=75)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_linewidth(10)
        ax.spines["left"].set_linewidth(10)

        # We change the fontsize of minor ticks label
        ax.tick_params(axis="both", which="major", labelsize=50)
        ax.tick_params(axis="both", which="minor", labelsize=50)

        if outdir:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            plt.savefig(
                "{}/{}_{}_alpha_{}.png".format(
                    outdir, sample_name, self.class_id, self.cluster_alpha
                ),
                bbox_inches="tight",
            )

        plt.show()


########################
# FUNCTION DEFINITIONS #
########################


def do_clustering(cell_type_df, eps_param, x_id_col, y_id_col, min_samples=1):
    """
    Performs DBSCAN clustering given a cell typing dataframe and returns list of labels of clusters that each cell belongs to.
    """

    n_cells = len(cell_type_df.index)
    print("n_cells = ", n_cells)

    # reset index so iteration on cell type df works:
    cell_type_df = cell_type_df.reset_index()
    cell_positions = cell_type_df[[x_id_col, y_id_col]].values.tolist()

    # perform DBSCAN clustering:
    clustering = DBSCAN(eps=eps_param, min_samples=min_samples).fit(cell_positions)

    return clustering.labels_


def df_select(df, cat, val):
    return df[df[cat] == val]


def process(
    data,
    image_id_col,
    x_id_col,
    y_id_col,
    clustering_id_col,
    cluster_class="Positive",
    cluster_eps=35,
    min_s=1,
):
    # read lobe data

    clustered_dfs = []

    for imagename in data[image_id_col].unique():
        # define empty dataframe for images with no clusters:
        no_cluster_df = pd.DataFrame(
            data={"Image": [imagename], "lesion_id": [-1], "cell_count": [0]}
        )

        image_data = df_select(data, image_id_col, imagename)

        if (
            len(image_data) > 1
        ):  # only proceed if there is more than one cell in the image
            positive_data = df_select(image_data, clustering_id_col, cluster_class)
            positive_data = positive_data.reset_index(drop=True)

            if (
                len(positive_data.index) > 0
            ):  # only proceed if there are positive cells in the image
                labels = do_clustering(
                    positive_data, cluster_eps, x_id_col, y_id_col, min_samples=min_s
                )
                positive_data.loc[:, "lesion_id"] = labels
                cluster_label, cluster_size = np.unique(labels, return_counts=True)
                cluster_id_df = pd.DataFrame(
                    data=zip(cluster_label, cluster_size),
                    columns=["Cluster ID", "Cluster N cells"],
                )
                merged = pd.merge(
                    positive_data,
                    cluster_id_df,
                    how="left",
                    left_on="lesion_id",
                    right_on="Cluster ID",
                )
                clustered_dfs.append(merged)
            else:
                print(f"Image {imagename} has no {cluster_class} cells, skipping...")
                clustered_dfs.append(no_cluster_df)

        else:
            print(f"Image {imagename} has only one cell, skipping...")
            clustered_dfs.append(no_cluster_df)

    clustered_data = pd.concat(clustered_dfs)

    return clustered_data


def measure_lesions(
    data,
    image_id_col,
    x_id_col,
    y_id_col,
    clustering_id_col,
    class_id,
    plot=False,
    outdir="./plots",
    sample_name="sample",
):
    """Function for measuring lesion size from a dataframe of cell coordinates.

    Args:
        data (pd.DataFrame):
            Dataframe of cell coordinates.
        image_id_col (str):
            Column name of image id, e.g. 'Image'.
        x_id_col (str):
            Column name of x coordinate, e.g. 'Centroid X µm'.
        y_id_col (str):
            Column name of y coordinate, e.g. 'Centroid Y µm'.
        clustering_id_col (str):
            Column name of cell class to find spatial clusters for, e.g. 'phenotype', 'class' etc.
        class_id (str):
            Class to find spatial clusters for, e.g. 'Positive', 'Negative' etc.
        plot (bool):
            Whether to plot the clusters. Default False. Produces plots with alphashapes for each cluster, the background tissue and points proportional to cluster sizes. This can be slow for large datasets.
        outdir (str):
            Output directory for plots. Default './plots'.

    Returns:
        summary (pd.DataFrame):
            Dataframe of lesion sizes, with columns 'Image', 'lesion_id', 'cell_count'.

            Example:
            Image	lesion_id	cell_count
            1179_22_1	0	1
            1179_22_1	1	56
            1179_22_1	2	8
            ...

    Example:
        >>> summary = measure_lesions(data, 'Image', 'Centroid X µm', 'Centroid Y µm', 'phenotype', 'Positive', plot = True, outdir = './plots')

    """

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    clustered_data = process(
        data,
        image_id_col=image_id_col,
        x_id_col=x_id_col,
        y_id_col=y_id_col,
        clustering_id_col=clustering_id_col,
        cluster_class=class_id,
    )

    if plot:
        plot_clusters_pub_proportional(
            clustered_data,
            data,
            x_id_col,
            y_id_col,
            sample_name,
            clustering_id_col,
            class_id,
            outdir,
            alphashape_param=0.05,
        )

    summary = (
        clustered_data.groupby([image_id_col, "lesion_id"])
        .size()
        .reset_index(name="cell_count")
    )

    points = data[[x_id_col, y_id_col]].values
    area = alphashape_area(points)

    summary["lobe_area"] = area
    summary["total_cell_area"] = data["Cell: Area"].sum()
    summary["total_lesions"] = len(summary)
    summary["lesions_per_unit_lobe_area"] = (
        summary["total_lesions"][0] / summary["lobe_area"][0]
    ) * 1e6
    summary["lesions_per_unit_cell_area"] = (
        summary["total_lesions"][0] / summary["total_cell_area"][0]
    ) * 1e6
    return summary


def alphashape_area(points, alphashape_param=0.03):
    """
    Function to calculate the area of a set of points with the alpha shape algorithm.
    """

    # create alphashape:
    print("Calculating alphashape...")
    alpha_shape = alphashape.alphashape(points, alphashape_param)
    area = alpha_shape.area
    print("Done.")
    return area


def main(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    flist = glob.glob(os.path.join(args.datapath, f"*{args.ext}"))

    all_data = []

    for f in flist:
        print(f"Clustering file: {f}...")

        fname = os.path.basename(f)

        data = pd.read_csv(f, sep="\t")

        clustered_data = process(
            data,
            image_id_col=args.image_id_col,
            x_id_col=args.x_coordname,
            y_id_col=args.y_coordname,
            clustering_id_col="Name",
            cluster_class="Positive",
        )

        if args.plot:
            plot_clusters_pub_proportional(
                clustered_data,
                data,
                args.x_coordname,
                args.y_coordname,
                fname.replace(args.ext, ""),
                "lesion_id",
                "Positive",
                f"{args.outdir}/plots",
                alphashape_param=0.05,
            )

        summary = (
            clustered_data.groupby(["Image", "lesion_id"])
            .size()
            .reset_index(name="cell_count")
        )
        summary.loc[summary["lesion_id"] == -1, ["cell_count"]] = np.nan

        summary.to_csv(
            os.path.join(args.outdir, fname.replace(args.ext, f"_cluster_sizes.csv")),
            sep="\t",
            index=False,
        )

        all_data.append(summary)

    # create summaries:
    all_data = pd.concat(all_data)
    all_data.to_csv(
        os.path.join(args.outdir, f"cluster_sizes_collated.csv"), sep="\t", index=False
    )

    all_data.fillna(0, inplace=True)
    all_data.groupby(["Image"])["cell_count"].agg(
        ["mean", "median", "std", "sum", "min", "max"]
    ).reset_index().to_csv(
        os.path.join(args.outdir, f"cluster_sizes_collated_summary.csv"),
        sep="\t",
        index=False,
    )

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Measure lesion sizes from cell coordinates."
    )
    parser.add_argument(
        "-d",
        "--datapath",
        type=str,
        help="Path to directory containing cell coordinate files.",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--ext",
        type=str,
        help="File extension of cell coordinate files.",
        default=".csv",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--image_id_col",
        type=str,
        help='Column name of image id, e.g. "Image".',
        default="Image",
        required=True,
    )
    parser.add_argument(
        "-x",
        "--x_coordname",
        type=str,
        help='Column name of x coordinate, e.g. "Centroid X µm".',
        default="Centroid X µm",
        required=True,
    )
    parser.add_argument(
        "-y",
        "--y_coordname",
        type=str,
        help='Column name of y coordinate, e.g. "Centroid Y µm".',
        default="Centroid Y µm",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Output directory for plots.",
        default="./results",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--plot",
        type=bool,
        help="Whether to plot the clusters. Default False.",
        default=False,
        required=False,
    )
    parser.add_argument(
        "-c",
        "--clustering_id_col",
        type=str,
        help='Column name of cell class to find spatial clusters for, e.g. "phenotype", "class" etc.',
        default="phenotype",
        required=True,
    )
    parser.add_argument(
        "--class",
        type=str,
        help='Class to find spatial clusters for, e.g. "Positive", "Negative" etc.',
        default="Positive",
        required=True,
    )
    args = parser.parse_args()
    main(args)
