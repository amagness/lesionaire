import seaborn as sns
import matplotlib.pyplot as plt
import os
import re
import matplotlib.pyplot as plt
import alphashape
from descartes import PolygonPatch
from shapely.geometry import MultiPolygon, Point, Polygon
import numpy as np


def plot_clusters_pub_proportional(clustered_df, 
                                allcells, 
                                xcoord_ref, 
                                ycoord_ref, 
                                sample_name,
                                lesion_id, 
                                clustering_cell_type, 
                                outdir, 
                                alphashape_param = 0.05):
    """
    Plots spatial clusters of cells with an alphashape pertaining to each spatial cluster. 
    Counts all other cell types which lie within alphashape and outputs dataframe.

    args:
        clustered_df: dataframe of clustered cells
        allcells: dataframe of all cells
        xcoord_ref: column name of x coordinates
        ycoord_ref: column name of y coordinates
        image_shape: shape of image

    """
    if os.path.exists(outdir) != True:
        os.makedirs(outdir)

    sns.set_style('white')
    alphashape_param = 0.03#0.04
    
    ## make alphashape and plot outline of lung
    # get points of cluster:
    X = allcells[xcoord_ref].values
    Y = allcells[ycoord_ref].values       
    points = list(zip(X, Y))
    lung_alpha_shape = alphashape.alphashape(points, alphashape_param)


    image_shape = (allcells[ycoord_ref].max(), allcells[xcoord_ref].max())
    
    # number of cells to split small and large clusters:
    cutoff = 5

    # create figure:
    scale = 5e-3
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (image_shape[1]*scale,image_shape[0]*scale))
    
    # plot lung:
    ax.add_patch(PolygonPatch(lung_alpha_shape, alpha=0.3)).set(facecolor='#f0b9e7', edgecolor='k')
        
    ## plot small clusters (<5 cells) as grey:
    small_clusters = clustered_df[clustered_df['Cluster N cells'] < cutoff]
    ax.scatter(small_clusters[xcoord_ref].values, small_clusters[ycoord_ref].values, alpha=0.1, color='k', s=small_clusters['Cluster N cells'].values*250)
    ax.scatter(small_clusters[xcoord_ref].values, small_clusters[ycoord_ref].values, alpha=1, s=1, color='k')

    large_clusters = clustered_df[clustered_df['Cluster N cells'] >= cutoff]

    # get unique cluster labels:
    unique_cluster_labels = large_clusters['Cluster ID'].unique()

    # loop through labels:
    for label in unique_cluster_labels:

        # get df for each unique label:
        cluster_df = large_clusters[large_clusters[lesion_id] == label]

        # get points of cluster:
        X = cluster_df[xcoord_ref].values
        Y = cluster_df[ycoord_ref].values       
        points = list(zip(X, Y))

        # only proceed if cells exist:
        if len(X) > 0:

            if label == -1:
                ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.75, color='k', s=50)

            else:
                # create alphashape:
                alpha_shape = alphashape.alphashape(points, alphashape_param)

                if alpha_shape.geom_type in ['Polygon', 'MultiPolygon']: # only process Polygon and Multipolygons i.e. ignore lines of cells which cannot contain other cells

                    # plot points and add patch of alphashape:
                    propscatter = ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, s=cluster_df['Cluster N cells'].values*250, alpha=0.05)
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=1, s=1, color='k')
                    ax.add_patch(PolygonPatch(alpha_shape, alpha=0.3)).set(facecolor='k', edgecolor='w')
#                     ax.legend()

                else:
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.1, color='k', s=150)
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=1, color='k', s=1)

    # update plot with title etc and save:
    title = '{}_{}_alpha_{}.png'.format(sample_name, clustering_cell_type, alphashape_param)
    ax.set_xlabel('Cell position X µm', fontsize=75)
    ax.set_ylabel('Cell position Y µm', fontsize=75)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(10)
    ax.spines['left'].set_linewidth(10)
    
    # We change the fontsize of minor ticks label 
    ax.tick_params(axis='both', which='major', labelsize=50)
    ax.tick_params(axis='both', which='minor', labelsize=50)

    #     plt.show()
    plt.savefig('{}/{}_{}_alpha_{}.png'.format(outdir,sample_name, clustering_cell_type, alphashape_param), bbox_inches = "tight")
    plt.close()


import re
import matplotlib.pyplot as plt
import alphashape
from descartes import PolygonPatch
from shapely.geometry import MultiPolygon, Point, Polygon

def assign_condition(mouseID):
    mouse_n = int(re.findall(r'\d+', mouseID)[0])
    if mouse_n <= 10:
        condition = 'PBS'
    elif (mouse_n >= 11) & (mouse_n < 20):
        condition = '5µg'
    else:
        condition = '50µg'
    return condition

def calc_lung_area(df):
    ''' Calculate lung area as the sum of all cell areas'''
    area = df['Cell: Area'].sum()
    return area

def df_select(df, cat, val):
    return df[df[cat]==val]

def plot_clusters(clustered_df, cluster_labels, xcoord_ref, ycoord_ref, image_shape, sample_name, clustering_cell_type, outdir, alphashape_param = 0.05):
    """
    Plots spatial clusters of cells with an alphashape pertaining to each spatial cluster. 
    Counts all other cell types which lie within alphashape and outputs dataframe.
    """
    if os.path.exists(outdir) != True:
        os.makedirs(outdir)

    # create figure:
    scale = 5e-3
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (image_shape[1]*scale,image_shape[0]*scale))

    # get unique cluster labels:
    unique_cluster_labels = np.unique(cluster_labels)

    # loop through labels:
    for label in unique_cluster_labels:
        
        
            

        # get df for each unique label:
        cluster_df = clustered_df[clustered_df['dbscan_cluster'] == label]

        # get points of cluster:
        X = cluster_df[xcoord_ref].values
        Y = cluster_df[ycoord_ref].values       
        points = list(zip(X, Y))

        # only proceed if cells exist:
        if len(X) > 0:
            
            if label == -1:
                ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.25, color='k', s=15)
            
            else:
                # create alphashape:
                alpha_shape = alphashape.alphashape(points, alphashape_param)

                if alpha_shape.geom_type in ['Polygon', 'MultiPolygon']: # only process Polygon and Multipolygons i.e. ignore lines of cells which cannot contain other cells

                    # plot points and add patch of alphashape:
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.5, s=20)
                    ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))

                else:
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.1, color='b', s=15)

    # update plot with title etc and save:
    title = '{}_{}_alpha_{}.png'.format(sample_name, clustering_cell_type, alphashape_param)
    ax.set_title(title, fontsize=18)
    ax.set_xlabel('Centroid X µm')
    ax.set_ylabel('Centroid Y µm')
    ax.set_ylim(0,image_shape[0])
    ax.set_xlim(0,image_shape[1])
#     plt.show()
    plt.savefig('{}/{}_{}_alpha_{}.png'.format(outdir,sample_name, clustering_cell_type, alphashape_param))
    plt.close()
    
def plot_clusters_pub(clustered_df, cluster_labels, xcoord_ref, ycoord_ref, image_shape, sample_name, clustering_cell_type, outdir, alphashape_param = 0.05):
    """
    Plots spatial clusters of cells with an alphashape pertaining to each spatial cluster. 
    Counts all other cell types which lie within alphashape and outputs dataframe.
    """
    if os.path.exists(outdir) != True:
        os.makedirs(outdir)

    # create figure:
    scale = 5e-3
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (image_shape[1]*scale,image_shape[0]*scale))

    # get unique cluster labels:
    unique_cluster_labels = np.unique(cluster_labels)

    # loop through labels:
    for label in unique_cluster_labels:
        
        
            

        # get df for each unique label:
        cluster_df = clustered_df[clustered_df['dbscan_cluster'] == label]

        # get points of cluster:
        X = cluster_df[xcoord_ref].values
        Y = cluster_df[ycoord_ref].values       
        points = list(zip(X, Y))

        # only proceed if cells exist:
        if len(X) > 0:
            
            if label == -1:
                ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.75, color='k', s=50)
            
            else:
                # create alphashape:
                alpha_shape = alphashape.alphashape(points, alphashape_param)

                if alpha_shape.geom_type in ['Polygon', 'MultiPolygon']: # only process Polygon and Multipolygons i.e. ignore lines of cells which cannot contain other cells

                    # plot points and add patch of alphashape:
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.5, s=200)
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=1, s=1, color='k')
                    ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))

                else:
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=0.1, color='k', s=150)
                    ax.scatter(cluster_df[xcoord_ref].values, cluster_df[ycoord_ref].values, alpha=1, color='k', s=1)

    # update plot with title etc and save:
    title = '{}_{}_alpha_{}.png'.format(sample_name, clustering_cell_type, alphashape_param)
    ax.set_title(title, fontsize=18)
    ax.set_xlabel('Centroid X µm')
    ax.set_ylabel('Centroid Y µm')
    ax.set_ylim(0,image_shape[0])
    ax.set_xlim(0,image_shape[1])
#     plt.show()
    plt.savefig('{}/{}_{}_alpha_{}.png'.format(outdir,sample_name, clustering_cell_type, alphashape_param))
    plt.close()