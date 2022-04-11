# genTools
Module for reading in genetic data, and processing it for dimensionality reduction methods, includes methods for reading several different formats of genetic data files and ancestry files.

## Usage
1. The following method can be used to read in a genotype data file and an ancestry call file, generating a masked matrix for the required ancestry:
    ```python
    get_masked_matrix(beagle_filename, vcf_filename, beagle_or_vcf, is_masked, vit_filename, fbk_filename, tsv_filename, vit_or_fbk_or_tsv, fb_or_msp, num_ancestries, ancestry, average_parents, prob_thresh)
    ```
    This is a wrapper method that calls the format-specific methods for reading files. 
    The arguments for this function are as follows:
    * `beagle_filename` (str): path to the Beagle file if the genotype data file is a Beagle file.
    * `vcf_filename` (str): path to the VCF file if the genotype data file is a VCF file.
    * `beagle_or_vcf` (int): `1` if the genotype data file is a Beagle file, or `2` if it is a VCF file.
    * `is_masked` (bool): `True` if an ancestry file is passed for ancestry-specific masking, or `False` otherwise.
    * `vit_filename` (str): path to the VIT file if the ancestry file is a VIT file.
    * `fbk_filename` (str): path to the FBK file if the ancestry file is a FBK file.
    * `tsv_filename` (str): path to the MSP.TSV or FB.TSV file if the ancestry file is a MSP.TSV or FB.TSV file.
    * `vit_or_fbk_or_tsv` (int): `1` if the ancestry file is a VIT file, `2` if it is an FBK file, or `3` if it is a MSP.TSV or FB.TSV file.
    * `fb_or_msp` (int): `1` if the TSV ancestry file is an FB.TSV file, or `2` if it is an MSP.TSV file.
    * `num_ancestries` (int): the total number of ancestries in the ancestry file.
    * `ancestry` (int): ancestry number of the ancestry for which dimensionality reduction is to be performed. Ancestry counter starts at 0 if the ancestry file is a MSP.TSV or FB.TSV file, and starts at 1 if it is a VIT or an FBK file.
    * `average_parents` (bool): `True` if the two parental haplotypes are to be combined (averaged) for each individual, or `False` otherwise.
    * `prob_thresh` (float): minimum probability threshold for a SNP to be assigned to an ancestry, (if the ancestry file is an FBK file or an FB.TSV file.)

    The function returns the following objects:
    * `masked_matrix` (matrix): the genotype matrix with the specified ancestry masked.
    * `ind_IDs` (vector): the vector of individual IDs.
    * `rs_IDs` (vector): the vector of rsIDs.

2. The following method can be used to read in a masked matrix and a labels file (with weights optional), and assign labels and weights to individuals or combine individuals, as specified:
    ```python
    process_labels_weights(labels_file, masked_matrix, ind_IDs, average_parents, is_weighted, save_masked_matrix, masked_matrix_filename)
    ```
    The arguments for this function are as follows:
    * `labels_file` (str): path to the labels file. It should be a TSV file where the first column has header `indID` and contains the individual IDs, and the second column has header `label` and contains the labels for all individuals. If `is_weighted` is specified as `True`, then the file must have a third column that has header `weight` and contains the weights for all individuals.
    NOTE: Individuals with positive weights are weighted accordingly. Individuals with zero weight are removed. Negative weights are used to combine individuals and replace them with a single average individual. Provide a weight of `-1` to the first set of individuals to be combined, `-2` to the second set of individuals to be combined, and so on. Each set of individuals that is to be combined must have the same label.
    * `masked_matrix` (matrix): the masked genotype matrix.
    * `ind_IDs` (vector): the vector of individual IDs.
    * `average_parents` (bool): `True` if the parental haplotypes are to be combined (averaged) for each individual, or `False` otherwise.
    * `is_weighted` (bool): `True` if weights are provided in the labels file, or `False` otherwise.  
    * `save_masked_matrix` (bool): `True` if the masked matrix is to be saved as a binary file, or `False` otherwise.
    * `masked_matrix_filename` (str): path to the masked matrix file. The masked matrix is saved in this file.

    The function returns the following objects:
    * `masked_matrix` (matrix): the genotype matrix with the specified ancestry masked.
    * `ind_IDs` (vector): the vector of individual IDs.
    * `labels` (vector): the vector of labels.
    * `weights` (vector): the vector of weights.

NOTE: There are 2 acceptable formats for SNP indices in the Beagle file: 
1. rsid: `rs` followed by the id (integer). For example, `rs12345`.
2. position: chromosome number (integer) followed by `_`, followed by the position (integer). For example, `10_12345`.
