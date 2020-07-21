## About DORGE
![image](https://github.com/biocq/DORGE/blob/master/DORGE_logo.svg)

The advent of the various large-scale genomics data especially The Cancer Genome Atlas (TCGA) facilitates the systematic characterization of cancer driver genes in Pan-Cancer analysis. To integrate available orthogonal datasets from diverse resources, an efficient approach that unbiasedly incorporates various types of features is still needed. Although other existing methods can also identify cancer driver genes, they cannot identify TSGs and OGs separately. The collection of features in previous approaches is typically limited and does not fully utilize genomic and epigenetic features that have been shown to effectively identify cancer driver genes in the past few years. To meet that need, we propose the tool DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes, a comprehensive machine-learning framework to discover cancer driver genes by integrating genetic and epigenetic data in a pan-cancer analysis.

Using the website, you can identify the cancer driver genes (TSGs and OGs) based on confident molecular profiling. The default parameter provides a convinent cancer driver gene list with TSG or OG role. Furthermore, one can also customize the prediction  by altering the parameters.

## Resources

* Complete codes (Rmd notebook file) and data to test the DORGE machine-learning model used in the DORGE paper can be found in DORGE_tool_reproduce.zip. 
For your convenience, illustrative codes of the DORGE machine-learning model used in the DORGE paper can be found at https://biocq.github.io/DORGE/DORGE.html. An online video that explains the code is available at https://www.youtube.com/watch?v=Pk8ZqoHK8zk.
* Codes to reproduce the results in DORGE paper can be found at https://github.com/biocq/DORGE_codes
* Codes to visualize the prediction by Shiny app can be found at https://github.com/biocq/DORGE_shiny
* Codes to process evaluation data can be found at the Evaluation folder.


## Citation

Please cite the paper (DORGE: Discovery of Oncogenes and Tumor SuppressoR Genes Using Genetic and Epigenetic Features, now submitted to Science Advances) if the resources are used elsewhere.

## Changelog
*  April 4. Update the README page.
*  April 27. Update the README page.
*  Jul 20. Update the README page.

## Issues

If you encounter any problems, please [file an issue](https://github.com/biocq/DORGE/issues) along with a detailed description.
