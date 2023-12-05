<div align="center" id="top"> 
  <img src="./docs/images/lesionaire.png" alt="Lesionaire" />

  &#xa0;

  <!-- <a href="https://lesionaire.netlify.app">Demo</a> -->
</div>

<h1 align="center">Lesionaire</h1>

<p align="center">
  <img alt="Github top language" src="https://img.shields.io/github/languages/top/amagness/lesionaire?color=56BEB8">

  <img alt="Github language count" src="https://img.shields.io/github/languages/count/amagness/lesionaire?color=56BEB8">

  <img alt="Repository size" src="https://img.shields.io/github/repo-size/amagness/lesionaire?color=56BEB8">

  <img alt="License" src="https://img.shields.io/github/license/amagness/lesionaire?color=56BEB8">

  <!-- <img alt="Github issues" src="https://img.shields.io/github/issues/amagness/lesionaire?color=56BEB8" /> -->

  <!-- <img alt="Github forks" src="https://img.shields.io/github/forks/amagness/lesionaire?color=56BEB8" /> -->

  <!-- <img alt="Github stars" src="https://img.shields.io/github/stars/amagness/lesionaire?color=56BEB8" /> -->
</p>

<!-- Status -->

<!-- <h4 align="center"> 
	🚧  Lesionaire 🚀 Under construction...  🚧
</h4> 

<hr> -->

<p align="center">
  <a href="#dart-about">About</a> &#xa0; | &#xa0; 
  <a href="#sparkles-features">Features</a> &#xa0; | &#xa0;
  <a href="#white_check_mark-requirements">Requirements</a> &#xa0; | &#xa0;
  <a href="#starting">Starting</a> &#xa0; | &#xa0;
  <a href="#license">License</a> &#xa0; | &#xa0;
  <a href="https://github.com/amagness" target="_blank">Author</a>
</p>

<br>

## About ##

Lesionaire is a small Python library for lesion size analysis and plotting with single cell spatial pathology data. The Lesionaire documentation is currently under development.

## Features ##

:heavy_check_mark: Easily quantify lesion size in cell coordinate data;\
:heavy_check_mark: Automatically determine boundaries of lesions;\
:heavy_check_mark: Produce informative summaries and plots;

## Requirements ##

Lesionaire requires 
```python
numpy 
pandas
scikit-learn, 
matplotlib 
seaborn
alphashape
descartes 
shapely
```

## Getting started ##

```bash
# Clone this project
$ git clone https://github.com/amagness/lesionaire

```

## Usage ##

```python
import lesionaire
import pandas as pd
import os

# read in the data
data = pd.read_csv('../data/test_data.txt', sep='\t')

# specify the column names for the x and y coordinates of the cells
x = 'Centroid X µm'
y = 'Centroid Y µm'
# specify the column names for the clustering id, the cluster class, the image id and the cell area. If you have use QuPath to identify positive and negative cells in a whole slide image, these will be 'Class' and 'Positive' respectively.
clustering_id_col = 'Class'
cluster_class = 'Positive'
image_id_col='Image' # if you have multiple images in your data file
cell_area_col = 'Cell: Area' # if you have cell area data (optional)

# create a lesionaire object
lesions = lesionaire.lesionData(data, 
                                x, 
                                y, 
                                clustering_id_col,, 
                                cluster_class, 
                                image_id_col = image_id_col, 
                                cell_area_col = cell_area_col)

imagename = os.path.split(f)[1].split('.')[0]
lesions.plot_lesions(outdir='../plots', sample_name=imagename)

```

## License ##

This project is under license from MIT. For more details, see the [LICENSE](LICENSE.md) file.


Lesionaire is developed by <a href="https://github.com/amagness" target="_blank">Alastair Magness</a> at the Francis Crick Institute.


