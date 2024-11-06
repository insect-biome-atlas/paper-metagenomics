# paper-metagenomics

This repo contains data and code used for matching metagenomic and metabarcoding results for the 40 IBA samples analyzed by Lopez-Clinton et al.

The code is in the R script `match_results.R`.

It relies on reading in metagenomics data from the file `metagenomics_genus_occurrence.tsv`, which is provided in the repo. We also provide the raw export from the metagenomics analysis in the file `metagenomics_binary_genus_results_onlybothdb_updated_12.10.2024.tsv`. The former file is based on the latter but includes additional taxonomic information from manual searches in NCBI taxonomy.

Furthermore, the `match_results.R` script relies on reading in the IBA metagenomic data using code in the [`utils` repo of the IBA organization on GitHub](https://github.com/insect-biome-atlas/utils/). This code needs to be pulled down from GitHub, and the link to it set correctly in the `match_results.R` script.

Finally, the script requires you to set the local paths to the IBA data, which you can download from FigShare following the links in the IBA data paper.

For more information on the IBA data, where to find it and how it is structured, see the [IBA data paper preprint](https://www.biorxiv.org/content/10.1101/2024.10.24.619818v1)

For the detailed analysis of the taxonomic annotations for the Lepidoptera clusters in the Swedish IBA data, see the [Iwaszkiewicz et al. preprint](https://www.biorxiv.org/content/10.1101/2024.10.25.620209v1)
