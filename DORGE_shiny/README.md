# DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes

<p>
	<a href="#">
	   <img src="https://travis-ci.org/jakelever/cancermine.svg?branch=master" />
	</a>
	<a href="#">
	   <img src="https://img.shields.io/badge/data-viewer-9e42f4.svg" />
	</a>
</p>

The advent of the various large-scale Genomics data especially The Cancer Genome Atlas (TCGA) facilitates the systematic characterization of cancer driver genes in Pan-Cancer analysis. To integrate available orthogonal datasets from diverse resources and different statistical backgrounds, an efficient approach that unbiasedly incorporates various types of features is still needed. Though other existing methods can also identify cancer driver genes but cannot identify TSGs and OGs separately. The collection of features in previous approaches is typically limited and does not fully utilize Genomic and Epigenetic features that have been evaluated in cancer driver genes in the past few years. To meet that need, we propose to develop Discovery of Oncogenes and Tumor SuppressoR Genes (DORGE), a comprehensive machine-learning framework to discover cancer driver genes on the basis of integrated genetic and epigenetic data in a Pan-Cancer analysis.

## System Requirements

This is a R shiny project which has been tested on MacOS and R 3.6.0 (64 bit only) but should work on other operating systems and R versions.

## Installation Guide

We provide a website of the Shiny app without installation http://216.127.179.19:3838/DORGE_shiny/ or http://45.12.109.148:3838/DORGE_shiny/ (be patient, somewhat slow to show).

You can also clone this repo using Git or download the [ZIP file](https://github.com/biocq/DORGE/archive/master.zip) of it.

```
git clone https://github.com/biocq/DORGE.git
```
* R software (https://www.r-project.org/) is needed.
* The codes in [app.R](https://github.com/biocq/DORGE/blob/master/app.R) are the Shiny codes used for the [web viewer].
* Please change the path of the folder that contains the files and set the path in app.R to the actual path, for example: setwd("/your_home_directory/DORGE/DORGE_shiny").
* If you find this website is helpful, you can use the code for your own projects. The list of dependencies is found at the top of the [app.R] file.

## Shiny app illustration

Using Shiny app, you can see/do:

* A full table with all genes in human with ranking and related information.
* Tables and barchart of Epigenetic features in different celllines or samples for a specific gene.
* A table for all genes with highest expression in different celllines or samples.
* Customize prediction of TSGs and OGs by redefining weights of different features (Plan to add in future).

## Citation

Please cite the paper (DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes Using Distinct Genetic and Epigenetic Features, https://www.biorxiv.org/content/10.1101/2020.07.21.213702v1) if the resources are used elsewhere.

## Changelog

* March 9. Update the name of the tool and fix a few bugs.
* April 29. Update the name of the tool and fix a few bugs.

## Issues

If you encounter any problems, please [file an issue](https://github.com/biocq/DORGE_shiny/issues) along with a detailed description.
