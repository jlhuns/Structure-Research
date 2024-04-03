# bio_465_structure_function
BIO 465 - Capstone Project


## Setup
- download and install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) if you don't already have it
- add conda to PATH (if necessary)
- create a new conda environement with the target by running `conda create --name <env>` (replace `<env>` with your desired environment name)
- activate the conda environment `conda activate <env>`
- install the necessary conda packages manually if needed
   - The following script should install the appropriate packages: `conda install python requests numpy pandas biopython networkx matplotlib && conda install -c etetoolkit ete3`


<!-- [post about conditional requirements files](https://stackoverflow.com/questions/29222269/is-there-a-way-to-have-a-conditional-requirements-txt-file-for-my-python-applica) -->


## To Run
### Main Pipeline
The main script for running the pipeline is found in `run_pipeline.py`. This file can be run using the commandline arguments of the filepath of a csv file containing a list of KEGG organism codes for target organisms and a KEGG id for a target KEGG Ortholog.
   - e.g. `python structure_function_pipeline.py "./target_organisms.csv" "K03841"` IS THIS IN A DIFFERENT ORDER?


#### Inputs
The `launch.json` file contains a configuration that should run the pipeline for KEGG Ortholog `K00937` with the organism list found in  `diverse_target_prokaryotes.csv`. This file can be updated to point to your desired KO ID and organism list, or you can simply run the pipeline from the command line.

   ##### KEGG Ortholog ID
   The first arguement that will be passed into the pipeline is the [KEGG Ortholog ID](https://www.genome.jp/kegg/ko.html) of the ortholog that you want to study. These ID's are formatted as 'K0#####'. This project focused on orthologs in metabolic pathways, as these are most likely to have conserved orthologs across a broad range of organisms.
   
   Example orthologs to use for a specific pathway can be found by:
   1. Visiting [KEGG Pathways](https://www.genome.jp/kegg/pathway.html) and selecting a pathway of choice (e.g. [map01200 - 'Carbon Metabolism'](https://www.genome.jp/pathway/map01200))
   2. Changing the URL to reflect a map specific to an appropriate organism. (e.g. 'https://www.genome.jp/pathway/map01200' --> https://www.genome.jp/pathway/eco01200 for E. Coli)
      - (Note that the pipeline will work best when the organism that you select is found in the list of target organisms that is provided by the sceond command-line arguement.)
   3. Selecting a green highlighted arrow or node, representing a protein (e.g. [glucokinase](https://www.genome.jp/entry/eco:b2388))
   4. The KO value for that protein is listed in the 'Name' box near the top of the information table (e.g. [K00845](https://www.genome.jp/entry/K00845) for glucokinase)

   ##### Organism List
   The second arguement passed into the pipeline is a CSV file containing information about which organisms should be considered when evaluating orthologs and conserved residues. The only requirement of the CSV file is a column titled `kegg_organism_code`, which contains the 3-letter code that KEGG uses to identify organisms (e.g. `eco` for E. Coli).

   The example `diverse_target_prokaryotes.csv` provided with the project is a list of 50 prokaryotes that have a significant amount of phylogenetic diversity.

### Pieces of the Pipeline

#### Ortholog Finder

#### Sequence Alignment

#### 3D-Clustering

### Figures

#### Figure 1

Make sure you have the ete3 python package installed

To replicate our Phylogenetic Tree Figures follow these steps:
1. Go to either `figures/figure_1/phylo_tree_prokaryotes.py` or `figures/figure_1/phylo_tree_eukaryotes.py`. Respectively they are the files to generate our trees using prokaryotes and eukaryotes
2. Run the `main` function at the bottom of the file.
3. Go to `figures/figure_1/figure_1_results`
4. In this directory you will see the tree you have created as either `eukaryote_tree.png` or `prokaryote_tree.png`

To generate your own Figure 1 Phylogenetic Tree follow these steps:
1. Run the pipeline with desired K0 id and organisms.
2. Go to `figures/figure_1/phylo_tree.py`
3. Run only the function `get_species_names_and_uniprot_ids_from_uniprot_entries()`
4. A file should have appeared in the `figure_1` directory called `included_organisms.txt`. This file shows each organism that had enough data to be used in the pipeline.
5. The first column will be the name of the uniprot file in `datafiles/uniprot_entries` corresponding to that organism. The second column will be separated from the first by just whitespace, and it has the species names.
6. The species name for the purpose of this code normally is two words, the genus and species. However, some names in the uniprot entries do not match convention. This file uses a regex string to identify what the species name should be. To know what species name you should keep track of, we have explained the regex code below in case the need to modify it arises
    + `OS\s+((?:(?!subsp\.)[^(.])+)?\s*`
    + The function `get_species_names_and_uniprot_ids_from_uniprot_entries()` in `figures/figure_1/phylo_tree.py` includes the above regex string which becomes a pattern that matches with the species name in each uniprot entry file so the name can be extracted.
    + `OS` this part finds the line in the uniprot file that starts with OS
    + `\s+` will match any number of whitespace characters
    + `((?:(?!subsp\.)[^(.])+)?` this is the first capturing group. Everything inside will be extracted as the species name. The final `?` means it can appear zero or one time
      + `(?:(?!subsp\.)[^(.])` this is the non-capturing group. It's denoted by `?:` and is used to group together sections of text
      + `(?!subsp\.)` This is the negative lookahead denoted by `?!`. Once the search reaches `subsp.` it will terminate.
      + `[^(.]` This matches everything that is not `(` or `.`.
      +  `+` means match as many characters as possible between one and infinity.
      + Altogether, this part matches as many characters as possible and will stop when it reaches `subsp.`, `(` or `.`.
    + `\s*` means match as many whitespaces as possible between zero and infinity. This is not part of the capturing group, but helps isolate it.
7. Go online to the NCBI Taxonomy Browser. https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
8. For every species in `included_organisms.txt` look up the taxonomy of that organism by typing the species name into the search bar 
9. Keep track of the taxonomies of every organism as you go. 
10. Once you have all taxonomies, decide how you want to group the organisms. Since every alignment and K0 id is different, the way you group your organisms will differ based one which had sufficient data to be included in the alignment. We grouped ours by class in prokaryote_tree.png. These groupings are for color coding purposes. 
11. Once you have decided how you want to group them, open `figures/figure_1/phylo_tree.py` again.
12. Go to the `set_styles_for_each_group()` function. You will be modifying this function. 
13. Follow the detailed instructions at the top of the function for how to write this code.

#### Figure 2


#### Figure 3
To generate Figure 3, run the file `figures/figure_3/make_figure_3.py` with 2 commandline arguements: 
1) the target UniProt ID of the protein which you want to analyze
2) the filepath to a CSV file that contains the cluster data generated by the pipeline. This cluster data can be found in `/datafiles/cluster_data/` directory.