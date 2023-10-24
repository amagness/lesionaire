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

## :dart: About ##

Lesionaire is a small Python library for lesion size analysis and plotting with single cell spatial pathology data. The Lesionaire documentation is currently under development.

## :sparkles: Features ##

:heavy_check_mark: Easily quantify lesion size in cell coordinate data;\
:heavy_check_mark: Automatically determine boundaries of lesions;\
:heavy_check_mark: Produce summaries and plots;

## :white_check_mark: Requirements ##

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
data = pd.read_csv('../data/test_data.txt', sep='\t')

L = lesionData(data = data,
                image_id_col = 'Image',
                x_id_col = 'Centroid X µm',
                y_id_col = 'Centroid Y µm',
                clustering_id_col = 'Class',
                class_id = 'Positive')

L.find_lesions()

lesion_fig = L.plot_lesions()
lesion_fig.savefig('./plots/lesion_fig.png')
```

## License ##

This project is under license from MIT. For more details, see the [LICENSE](LICENSE.md) file.


Made with :heart: by <a href="https://github.com/amagness" target="_blank">Alastair Magness</a> at the Francis Crick Institute.

&#xa0;

<a href="#top">Back to top</a>
