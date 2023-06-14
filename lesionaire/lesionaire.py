import sys

sys.path.append('./')
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd
from lesionaire.plots import plot_clusters_pub_proportional

from sklearn.cluster import DBSCAN
from tqdm import *

sys.path.append("../")



########################
# FUNCTION DEFINITIONS #
########################

def do_clustering(cell_type_df, eps_param, xcoord_ref, ycoord_ref,min_samples=1):
    """
    Performs DBSCAN clustering given a cell typing dataframe and returns list of labels of clusters that each cell belongs to.
    """

    n_cells = len(cell_type_df.index)
    print('n_cells = ', n_cells)

    # reset index so iteration on cell type df works:
    cell_type_df = cell_type_df.reset_index()

    # create list of X,Y positions in [[X0,Y0], [X1,Y1]...[Xi,Yi]] format for clustering:
    cell_positions = []
    for i in range(n_cells):
        cell_positions.append([cell_type_df[xcoord_ref][i], cell_type_df[ycoord_ref][i]])

    # perform DBSCAN clustering:
    clustering = DBSCAN(eps=eps_param, min_samples=min_samples).fit(cell_positions)
    labels = clustering.labels_
    return labels

def df_select(df, cat, val):
    return df[df[cat]==val]

def process(data, image_id_col, x_id_col, y_id_col, clustering_id_col, 
            cluster_class = 'Positive', cluster_eps = 35, whole_lung_eps = 1000, min_s = 1):
    
    # read lobe data

    clustered_dfs = []

    

    for imagename in data[image_id_col].unique():
        
        # define empty dataframe for images with no clusters:
        no_cluster_df = pd.DataFrame(data={'Image': [imagename],'lesion_id':[-1],'cell_count':[0]})

        image_data = df_select(data, image_id_col, imagename)

        if len(image_data) > 1: # only proceed if there is more than one cell in the image
            positive_data = df_select(image_data,  clustering_id_col, cluster_class)
            positive_data = positive_data.reset_index(drop=True)

            if len(positive_data.index) > 0: # only proceed if there are positive cells in the image
                labels = do_clustering(positive_data, cluster_eps, x_id_col, y_id_col, min_samples=min_s)
                positive_data.loc[:,'lesion_id'] = labels
                cluster_label, cluster_size = np.unique(labels, return_counts=True)
                cluster_id_df = pd.DataFrame(data=zip(cluster_label, cluster_size), columns = ["Cluster ID", "Cluster N cells"])
                merged = pd.merge(positive_data, cluster_id_df, how='left', left_on='lesion_id', right_on='Cluster ID')
                clustered_dfs.append(merged)
            else:
                print(f'Image {imagename} has no {cluster_class} cells, skipping...')
                clustered_dfs.append(no_cluster_df)

        else:
            print(f'Image {imagename} has only one cell, skipping...')
            clustered_dfs.append(no_cluster_df)

    clustered_data = pd.concat(clustered_dfs)

    return clustered_data

class lesions:
    """_summary_
    """

    def __init__(self,data, x_id_col, y_id_col, clustering_id_col, class_id, image_id_col=None, lobe_id_col = None) -> None:

        """_summary_
        """
        
        self.data = data
        self.image_id_col = image_id_col
        self.lobe_id_col = lobe_id_col
        self.x_id_col = x_id_col
        self.y_id_col = y_id_col
        self.clustering_id_col = clustering_id_col
        self.class_id = class_id
        self.cluster_eps = 35
        self.cluster_alpha = 0.05
        self.lung_alpha = 0.05
        self.min_s = 1

    def find_lesions(self):
        """Find lesions with density-based clustering (DBSCAN).

        Returns:
            pd.DataFrame: Original data with lesion ID column added.
        """
        self.data = process(self.data, self.image_id_col, self.x_id_col, self.y_id_col, self.clustering_id_col, 
            cluster_class = self.class_id, cluster_eps = self.cluster_eps, min_s = self.min_s)
        return self.data
    
    def cell_counts(self):
        if self.lobe_id_col:
            self.counts = self.data.groupby([self.image_id_col,self.lobe_id_col,'lesion_id']).size().reset_index(name='cell_count')
        else:
            self.counts = self.data.groupby([self.image_id_col,'lesion_id']).size().reset_index(name='cell_count')
        return self.counts
    
    def summary(self):
        self.summary = self.cell_counts().groupby(self.image_id_col).agg({'cell_count':['mean','std','count']})
        self.summary.columns = ['_'.join(col).strip() for col in self.summary.columns.values]
        self.summary = self.summary.reset_index()
        return self.summary

def measure_lesions(data, image_id_col, x_id_col, y_id_col, clustering_id_col, class_id, plot = False, outdir = './plots'):

    ''' Function for measuring lesion size from a dataframe of cell coordinates.

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

    '''

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    clustered_data = process(data, image_id_col=image_id_col, x_id_col=x_id_col, 
                y_id_col=y_id_col, clustering_id_col=clustering_id_col, 
                cluster_class=class_id)
    
    if plot:
        plot_clusters_pub_proportional(clustered_data, 
                                data, 
                                x_id_col, 
                                y_id_col, 
                                clustering_id_col,
                                class_id, 
                                outdir, 
                                alphashape_param = 0.05)
        
    summary = clustered_data.groupby([image_id_col, 'lesion_id']).size().reset_index(name='cell_count')
    return summary

def main(args):



    datapath = '/camp/lab/swantonc/inputs/histopath/WH/1179_22/MZ analysis 2023-05-19/results'
    ext = '.tsv'
    image_id_col = 'Image'
    x_coordname = 'Centroid X µm'
    y_coordname = 'Centroid Y µm'
    outdir = '../results/1179_22/MZ analysis 2023-05-19'
    plot = True

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    flist = glob.glob(os.path.join(datapath, f'*{ext}'))

    all_data = []

    for f in flist:

        print(f"Clustering file: {f}...")

        fname = os.path.basename(f)

        data = pd.read_csv(f,sep='\t') 

        clustered_data = process(data, image_id_col=image_id_col, x_id_col=x_coordname, 
                y_id_col=y_coordname, clustering_id_col='Name', 
                cluster_class='Positive')
        
        
        if plot:
            plot_clusters_pub_proportional(clustered_data, 
                                data, 
                                x_coordname, 
                                y_coordname, 
                                fname.replace(ext, ''),
                                'lesion_id',
                                'Positive', 
                                f'{outdir}/plots', 
                                alphashape_param = 0.05)
        
        summary = clustered_data.groupby(['Image', 'lesion_id']).size().reset_index(name='cell_count')
        summary.loc[summary['lesion_id'] == -1, ['cell_count']] = np.nan

        summary.to_csv(os.path.join(outdir, fname.replace(ext, f'_cluster_sizes.csv')), sep='\t', index=False)

        all_data.append(summary)


    # create summaries:    
    all_data = pd.concat(all_data)
    all_data.to_csv(os.path.join(outdir, f'cluster_sizes_collated.csv'), sep='\t', index=False)

    all_data.fillna(0, inplace=True)
    all_data.groupby(['Image'])['cell_count'].agg(['mean', 'median', 'std', 'sum', 'min', 'max']).reset_index().to_csv(os.path.join(outdir, f'cluster_sizes_collated_summary.csv'), sep='\t', index=False)

    print('Done.')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Measure lesion sizes from cell coordinates.')
    parser.add_argument('-d', '--datapath', type=str, help='Path to directory containing cell coordinate files.', required=True)
    parser.add_argument('-e', '--ext', type=str, help='File extension of cell coordinate files.', default='.csv', required=True)
    parser.add_argument('-i', '--image_id_col', type=str, help='Column name of image id, e.g. "Image".', default='Image', required=True)
    parser.add_argument('-x', '--x_coordname', type=str, help='Column name of x coordinate, e.g. "Centroid X µm".', default='Centroid X µm', required=True)
    parser.add_argument('-y', '--y_coordname', type=str, help='Column name of y coordinate, e.g. "Centroid Y µm".', default='Centroid Y µm', required=True)
    parser.add_argument('-o', '--outdir', type=str, help='Output directory for plots.', default='./results', required=True)
    parser.add_argument('-p', '--plot', type=bool, help='Whether to plot the clusters. Default False.', default=False, required=False)
    parser.add_argument('-c', '--clustering_id_col', type=str, help='Column name of cell class to find spatial clusters for, e.g. "phenotype", "class" etc.', default='phenotype', required=True)
    parser.add_argument('--class', type=str, help='Class to find spatial clusters for, e.g. "Positive", "Negative" etc.', default='Positive', required=True)
    args = parser.parse_args()
    main(args)

            