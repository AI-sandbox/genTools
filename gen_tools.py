# Library Importation
import allel
import numpy as np
import pandas as pd
import sys
import time
import warnings
import logging
import datetime
from os import path
warnings.simplefilter("ignore", category=RuntimeWarning)


def process_vit(vit_file):
    """                                                                                       
    Viterbi File Processing                                                 
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    vit_file    : string with the path of our viterbi file.                                                                       

    Returns                                                                                   
    -------   
    ancestry_matrix  : (m, n) array
                  Viterbi Matrix indicating the ancestry for individual n at
                  position m.                                                                                                                                             
    """
    start_time = time.time()
    vit_matrix = []
    with open(vit_file) as file:
        for x in file:
            x_split = x.replace('\n', '').split('\t')
            vit_matrix.append(np.array(x_split[1:-1]))
    ancestry_matrix = np.stack(vit_matrix, axis=0).T
    logging.info("VIT Processing Time: --- %s seconds ---" % (time.time() - start_time))
    return ancestry_matrix

def process_fbk(fbk_file, num_ancestries, prob_thresh):    
    """                                                                                       
    FBK File Processing                                                 
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    fbk_file       : string with the path of our fbk file.    
    num_ancestries : number of distinct ancestries in dataset.
    prob_thres     : probability threshold for ancestry assignment.

    Returns                                                                                   
    -------   
    ancestry_matrix  : (m, n) array
                  Viterbi Matrix indicating the ancestry for individual n at
                  position m.                                                                                                                                             
    """
    start_time = time.time()
    df_fbk = pd.read_csv(fbk_file, sep=" ", header=None)
    fbk_matrix = df_fbk.values[:, :-1]
    ancestry_matrix = np.zeros((fbk_matrix.shape[0], int(fbk_matrix.shape[1] / num_ancestries)), dtype=np.int8)
    for i in range(num_ancestries):
        ancestry = i+1
        ancestry_matrix += (fbk_matrix[:, i::num_ancestries] > prob_thresh) * 1 * ancestry
    ancestry_matrix = ancestry_matrix.astype(str)
    logging.info("FBK Processing Time: --- %s seconds ---" % (time.time() - start_time))
    return ancestry_matrix

def process_tsv_fb(tsv_file, num_ancestries, prob_thresh, positions, gt_matrix, rs_IDs):
    """                                                                                       
    tsv_fb File Processing                                                 
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    tsv_file       : string with the path of our tsv file.    
    num_ancestries : number of distinct ancestries in dataset.
    prob_thres     : probability threshold for ancestry assignment.
    positions      :
    gt_matrix      : (m,n) array
                     Genome matrix indicating ref/alt letter for individual n
                     at position m.

    Returns                                                                                   
    -------   
    ancestry_matrix  : (m, n) array
                  Viterbi Matrix indicating the ancestry for individual n at
                  position m.                                                                                                                               
    """
    start_time = time.time()
    df_tsv = pd.read_csv(tsv_file, sep="\t", skiprows=1)
    tsv_positions = df_tsv['physical_position'].tolist()
    df_tsv.drop(columns = ['physical_position', 'chromosome', 'genetic_position', 'genetic_marker_index'], inplace=True)
    tsv_matrix = df_tsv.values
    i_start = positions.index(tsv_positions[0])
    if tsv_positions[-1] in positions:
        i_end = positions.index(tsv_positions[-1]) + 1
    else:
        i_end = len(positions)
    gt_matrix = gt_matrix[i_start:i_end, :]
    positions = positions[i_start:i_end]
    rs_IDs = rs_IDs[i_start:i_end]
    prob_matrix = np.zeros((len(positions), tsv_matrix.shape[1]), dtype=np.float16)

    i_tsv = -1
    next_pos_tsv = tsv_positions[i_tsv+1]
    for i in range(len(positions)):
        pos = positions[i]
        if pos >= next_pos_tsv and i_tsv + 1 < tsv_matrix.shape[0]:
            i_tsv += 1
            probs = tsv_matrix[i_tsv, :]
            if i_tsv + 1 < tsv_matrix.shape[0]:
                next_pos_tsv = tsv_positions[i_tsv+1]
        prob_matrix[i, :] = probs

    tsv_matrix = prob_matrix
    ancestry_matrix = np.zeros((tsv_matrix.shape[0], int(tsv_matrix.shape[1] / num_ancestries)), dtype=np.int8)
    for i in range(num_ancestries):
        ancestry = i+1
        ancestry_matrix += (tsv_matrix[:, i::num_ancestries] > prob_thresh) * 1 * ancestry
    ancestry_matrix -= 1
    ancestry_matrix = ancestry_matrix.astype(str)
    logging.info("TSV Processing Time: --- %s seconds ---" % (time.time() - start_time))
    return ancestry_matrix, gt_matrix, rs_IDs

def process_tsv_msp(tsv_file, positions, gt_matrix, rs_IDs):
    """                                                                                       
    tsv_msp File Processing                                                 
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    tsv_file       : string with the path of our tsv file.    
    positions      :
    gt_matrix      : (m,n) array
                     Genome matrix indicating ref/alt letter for individual n
                     at position m.

    Returns                                                                                   
    -------   
    ancestry_matrix  : (m, n) array
                  Viterbi Matrix indicating the ancestry for individual n at
                  position m.                                                                                                                               
    """
    start_time = time.time()
    df_tsv = pd.read_csv(tsv_file, sep="\t", skiprows=1)
    tsv_spos = df_tsv['spos'].tolist()
    tsv_epos = df_tsv['epos'].tolist()
    df_tsv.drop(columns = ['#chm', 'spos', 'epos', 'sgpos', 'egpos', 'n snps'], inplace=True)
    tsv_matrix = df_tsv.values
    i_start = positions.index(tsv_spos[0])
    if tsv_epos[-1] in positions:
        i_end = positions.index(tsv_epos[-1])
    else:
        i_end = len(positions)
    gt_matrix = gt_matrix[i_start:i_end, :]
    positions = positions[i_start:i_end]
    rs_IDs = rs_IDs[i_start:i_end]
    ancestry_matrix = np.zeros((len(positions), tsv_matrix.shape[1]), dtype=np.int8)

    i_tsv = -1
    next_pos_tsv = tsv_spos[i_tsv+1]
    for i in range(len(positions)):
        pos = positions[i]
        if pos >= next_pos_tsv and i_tsv + 1 < tsv_matrix.shape[0]:
            i_tsv += 1
            ancs = tsv_matrix[i_tsv, :]
            if i_tsv + 1 < tsv_matrix.shape[0]:
                next_pos_tsv = tsv_spos[i_tsv+1]
        ancestry_matrix[i, :] = ancs

    ancestry_matrix = ancestry_matrix.astype(str)
    logging.info("TSV Processing Time: --- %s seconds ---" % (time.time() - start_time))
    return ancestry_matrix, gt_matrix, rs_IDs

def process_beagle(beagle_file, rs_ID_dict, rsid_or_chrompos):
    """                                                                                       
    Beagle File Processing                                                 
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    beagle_file : string with the path of our beagle file.  
    rs_ID_dict  : dictionary showing the previous encoding for a specific 
                  rs ID.

    Returns                                                                                   
    -------   
    gt_matrix.  : (m, n) array
                  Genetic matrix indicating the encoding for individual n at 
                  poisition m. 
    ind_IDs     : (n,) array
                  Individual IDs for all individuals in the matrix. 
    rs_IDs      : (m,) array
                  rs IDs of all the positions included in our matrix. 
    rs_ID_dict  :
                  Encoding dictionary for each of the positions in dataset.                                                                                                                                               
    """
    start_time = time.time()
    rs_IDs = []
    lis_beagle = []
    with open(beagle_file) as file:
        x = file.readline()
        x_split = x.replace('\n', '').split('\t')
        ind_IDs = x_split[2:]
        ind_IDs = np.array(ind_IDs)
        for x in file:
            x_split = x.replace('\n', '').split('\t')
            if rsid_or_chrompos == 1:
                rs_IDs.append(int(x_split[1][2:]))
            elif rsid_or_chrompos == 2:
                rs_ID_split = x_split[1].split('_')
                rs_IDs.append(np.float64(rs_ID_split[0] + '.' + rs_ID_split[1][::-1]))
            else:
                sys.exit("Illegal value for rsid_or_chrompos. Choose 1 for rsID format or 2 for Chromosome_position format.")
            lis_beagle.append(x_split[2:])
    
    gt_matrix = np.zeros((len(lis_beagle),len(lis_beagle[0])), dtype=np.float16)
    
    processed_IDs = rs_ID_dict.keys()
    for i in range(len(lis_beagle)):
        # Check how we usually encode:
        if (rs_IDs[i] in processed_IDs):
            ref = rs_ID_dict[rs_IDs[i]]
        else:
            ref = lis_beagle[i][0]
            rs_ID_dict[rs_IDs[i]] = ref

        for j in range(1, len(lis_beagle[i])):
            gt_matrix[i, j] = (lis_beagle[i][j] != ref)*1

    logging.info("Beagle Processing Time: --- %s seconds ---" % (time.time() - start_time))

    return gt_matrix, ind_IDs, rs_IDs, rs_ID_dict

def process_vcf(vcf_file, rs_ID_dict, rsid_or_chrompos):
    """                                                                                       
    VCF File Processing                                                 
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    vcf_file    : string with the path of our vcf file. 
    rs_ID_dict  : dictionary showing the previous encoding for a specific 
                  rs ID.

    Returns                                                                                   
    -------   
    gt_matrix.  : (m, n) array
                  Genetic matrix indicating the encoding for individual n at 
                  poisition m. 
    ind_IDs     : (n,) array
                  Individual IDs for all individuals in the matrix. 
    rs_IDs      : (m,) array
                  rs IDs of all the positions included in our matrix. 
    positions   :  
    rs_ID_dict  : dictionary showing the previous encoding for a specific 
                  rs ID.
    """
    start_time = time.time()
    vcf = allel.read_vcf(vcf_file)
    gt = vcf['calldata/GT']
    n_variants, n_samples, ploidy = gt.shape
    gt_matrix = gt.reshape(n_variants, n_samples * ploidy).astype(np.float16)
    np.place(gt_matrix, gt_matrix < 0, np.nan)
    if rsid_or_chrompos == 1:
        IDs = vcf['variants/ID']
        rs_IDs = [int(x[2:]) for x in IDs]
    elif rsid_or_chrompos == 2:
        rs_IDs = []
        for i in range(len(vcf['variants/CHROM'])):
            rs_IDs.append(np.float64(vcf['variants/CHROM'][i] + '.' + str(vcf['variants/POS'][i])[::-1]))
    else:
        sys.exit("Illegal value for rsid_or_chrompos. Choose 1 for rsID format or 2 for Chromosome_position format.")
    ref_vcf = vcf['variants/REF']
    samples = vcf['samples']
    ind_IDs = []
    for sample in samples:
        ind_IDs.append(sample + '_A')
        ind_IDs.append(sample + '_B')
    ind_IDs = np.array(ind_IDs)
    positions = vcf['variants/POS'].tolist()
    
    processed_IDs = rs_ID_dict.keys()
    for i in range(len(rs_IDs)):
        rs_ID = rs_IDs[i]
        if (rs_ID in processed_IDs):
            ref = rs_ID_dict[rs_ID]
        else:
            ref = ref_vcf[i]
            rs_ID_dict[rs_ID] = ref
        if ref != ref_vcf[i]:
            gt_matrix[i, :] = 1 - gt_matrix[i, :]
    
    logging.info("VCF Processing Time: --- %s seconds ---" % (time.time() - start_time))
    return gt_matrix, ind_IDs, rs_IDs, positions, rs_ID_dict


###########################################################################
def mask(ancestry_matrix, gt_matrix, unique_ancestries, dict_ancestries, average_parents = False):
    """                                                                                       
    Masking Function for each of the available ancestries.                                               
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    ancestry_matrix   : (m, n) array
                        Ancestry Matrix indicating the ancestry for individual n at
                        position m. 
    gt_matrix         : (m, n) array
                        Genetic matrix indicating the encoding for individual n at 
                        poisition m. 
    unique_ancestries : list of u distinct unique ancestries in our ancestry file.
    average_parents   : Boolean to combine haplotypes for each individuals.

    Returns                                                                                   
    -------   
    masked_matrices : (m, n, u) dictionary/3D array
                      masked matrices for each of the distinct ancestries in the 
                      dataset.            
    """
    
    start_time = time.time()
    masked_matrices = {}
    for i in range(len(unique_ancestries)):
        ancestry = unique_ancestries[i]
        dict_ancestry = dict_ancestries[i]
        masked = np.empty(ancestry_matrix.shape[0] * ancestry_matrix.shape[1], dtype = np.float16)
        masked[:] = np.NaN
        arg = ancestry_matrix.reshape(-1) == ancestry
        masked[arg] = gt_matrix.reshape(-1)[arg]
        logging.info("Masking for ancestry " + str(ancestry) + " --- %s seconds ---" % (time.time() - start_time))

        if (average_parents == True):
            masked_matrices[dict_ancestry] = average_parent_snps(masked.reshape(ancestry_matrix.shape).astype(np.float16))    
        else:
            masked_matrices[dict_ancestry] = masked.reshape(ancestry_matrix.shape).astype(np.float16)
        start_time = time.time()
    return masked_matrices

def average_parent_snps(matrix):
    """                                                                                       
    Combining Haplotypes Function.                                               
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    matrix      : (m, n) array 
                  The masked matrix for an ancestry. 

    Returns                                                                                   
    -------   
    new_matrix  : (m, n/2) array 
                  The combined masked matrix with the actual number of individuals.
                                  
    """
    start = time.time()
    new_matrix = np.zeros((matrix.shape[0],int(matrix.shape[1]/2)), dtype = np.float16)
    for i in range(0,matrix.shape[1],2):
        new_matrix[:, int(i/2)] = np.nanmean(matrix[:,i:i+2],axis=1, dtype = np.float16)
    logging.info("Combining time --- %s seconds ---" % (time.time() - start))
    return new_matrix

###########################################################################
def get_masked_matrix(beagle_vcf_filename, beagle_or_vcf, vit_fbk_fbtsv_msptsv_filename, vit_or_fbk_or_fbtsv_or_msptsv,
                      is_mixed, is_masked, num_ancestries, average_parents, prob_thresh, rs_ID_dict, rsid_or_chrompos):
    """                                                                                       
    Input Parameter Parser                                                
                                                                                               
    Parameters                                                                                
    ----------       
    beagle_vcf_filename : string with the path of our beagle/vcf file.  
    beagle_or_vcf       : int
                          indicates file type 1=beagle, 2=vcf
                          
    vit_fbk_tsv_filename: string with the path of our vit/fbk/tsv file.
    vit_or_fbk_or_tsv   : int
                          indicates file type 1=vit, 2=fbk, 3=tsv
    fb_or_msp           : int
                          indicates 1=fb, 2=msp
    is_masked           : boolean
                          indicates if output matrix needs to be masked.
    num_ancestries      : int
                          number of distinct ancestries in dataset
    average_parents     : boolean
                          indicates whether to combine haplotypes for each individuals.
    prob_thresh         : float  
                          probability threshold for ancestry assignment.
    rs_ID_dict          : dictionary showing the previous encoding for a specific 
                          rs ID.

    Returns                                                                                   
    -------   
    masked_matrices  : (m, n)/(m, n, num_ancestries) array/dictionary
                       unmasked matrix/masked matrices for each of the distinct ancestries in the 
                       dataset.
    ind_IDs          : (n,) array
                       Individual IDs for all individuals in the matrix. 
    rs_IDs           : (m,) array
                       rs IDs of all the positions included in our matrix. 
    rs_ID_dict       :
                       Encoding dictionary for each of the positions in dataset.             
    """
    
    if beagle_or_vcf == 1:
        gt_matrix, ind_IDs, rs_IDs, rs_ID_dict = process_beagle(beagle_vcf_filename, rs_ID_dict, rsid_or_chrompos)
    elif beagle_or_vcf == 2:
        gt_matrix, ind_IDs, rs_IDs, positions, rs_ID_dict = process_vcf(beagle_vcf_filename, rs_ID_dict, rsid_or_chrompos)
        
    if is_masked and vit_or_fbk_or_fbtsv_or_msptsv != 0:
        if vit_or_fbk_or_fbtsv_or_msptsv == 1:
            ancestry_matrix = process_vit(vit_fbk_fbtsv_msptsv_filename)
        elif vit_or_fbk_or_fbtsv_or_msptsv == 2:
            ancestry_matrix = process_fbk(vit_fbk_fbtsv_msptsv_filename, num_ancestries, prob_thresh)
        elif vit_or_fbk_or_fbtsv_or_msptsv == 3:
            ancestry_matrix, gt_matrix, rs_IDs = process_tsv_fb(vit_fbk_fbtsv_msptsv_filename, num_ancestries, prob_thresh,
                                                                positions, gt_matrix, rs_IDs)
        elif vit_or_fbk_or_fbtsv_or_msptsv == 4:
            ancestry_matrix, gt_matrix, rs_IDs = process_tsv_msp(vit_fbk_fbtsv_msptsv_filename, positions, gt_matrix, rs_IDs)

        if vit_or_fbk_or_fbtsv_or_msptsv == 1 or vit_or_fbk_or_fbtsv_or_msptsv == 2:
            unique_ancestries = [str(i) for i in np.arange(1, num_ancestries+1)]
        else:
            unique_ancestries = [str(i) for i in np.arange(0, num_ancestries)]
        if is_mixed:
            dict_ancestries = [str(i) for i in np.arange(0, num_ancestries)]
        else:
            dict_ancestries = unique_ancestries
        masked_matrices = mask(ancestry_matrix, gt_matrix, unique_ancestries, dict_ancestries, average_parents)
    
    else:
        if not is_masked:
            dict_ancestries = [str(i) for i in np.arange(1, num_ancestries+1)]
        elif is_mixed or beagle_or_vcf == 2:
            dict_ancestries = [str(i) for i in np.arange(0, num_ancestries)]
        else:
            dict_ancestries = [str(i) for i in np.arange(1, num_ancestries+1)]
        masked_matrices = {}
        if average_parents:
            gt_matrix_avg = average_parent_snps(gt_matrix)
            for ancestry in dict_ancestries:
                masked_matrices[ancestry] = gt_matrix_avg
        else:
            for ancestry in dict_ancestries:
                masked_matrices[ancestry] = gt_matrix
        logging.info("No masking")
        
    return masked_matrices, ind_IDs, rs_IDs, rs_ID_dict


def array_process(root_dir, beagle_vcf_file, vit_fbk_fbtsv_msptsv_file, num_arrays, num_ancestries, average_parents, prob_thresh, is_masked, rsid_or_chrompos): 
    """                                                                                       
    Dataset processing of each of the individual arrays.                                               
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    root_dir         : string                                                                         
                       Directory of array folders containing input files.
    beagle_vcf_file  : string
                       Beagle/VCF Filename defined by user.
    vit_fbk_tsv_file : string
                       Viterbi/TSV/FBK Filename defined by user.
    beagle_or_vcf    : int
                       indicates file type 1=beagle, 2=vcf
    vit_or_fbk_or_tsv: int
                       indicates file type 1=vit, 2=fbk, 3=tsv
    fb_or_msp        : int
                       indicates 1=fb, 2=msp     
    num_arrays       : Total number of arrays in dataset.
    num_ancestries   : Number of unique ancestries in dataset. 
    average_parents  : boolean
                       Indicates whether to combine haplotypes for each individual.
    prob_thresh      : float  
                       Probability threshold for ancestry assignment.
    is_masked        : boolean
                       indicates if output matrix needs to be masked. 
    save_masks       : boolean
                       indicates if mask files needs to be saved.
    masks_file       : string
                       npz filename defined by user to save the mask files.
                                                                                   
    Returns                                                                                   
    -------                                                                                   
    masks      : (num_arrays, ) list                                                                          
                 List of masked matrices for each ancestries at each given array.
    rs_ID_list : (num_arrays, ) list 
                 List of rs IDs for each of the processed arrays.
    ind_ID_list: 
                 List of individual IDs for each of the processed arrays.

    """
    beagle_or_vcf_list = []
    vit_or_fbk_or_fbtsv_or_msptsv_list = []
    beagle_vcf_file_list = []
    vit_fbk_fbtsv_msptsv_file_list = []
    for i in range(1, num_arrays+1):
        beagle_vcf_file_path = root_dir + "array" + str(i) + '/' + beagle_vcf_file
        if path.exists(beagle_vcf_file_path + '.beagle'):
            beagle_vcf_file_list.append(beagle_vcf_file_path + '.beagle')
            beagle_or_vcf_list.append(1)
        elif path.exists(beagle_vcf_file_path + '.vcf'):
            beagle_vcf_file_list.append(beagle_vcf_file_path + '.vcf')
            beagle_or_vcf_list.append(2)
        else:
            sys.exit("No beagle/vcf file exists with this name in array " + str(i))
        vit_fbk_fbtsv_msptsv_file_path = root_dir + "array" + str(i) + '/' + vit_fbk_fbtsv_msptsv_file
        if path.exists(vit_fbk_fbtsv_msptsv_file_path + '.vit'):
            vit_fbk_fbtsv_msptsv_file_list.append(vit_fbk_fbtsv_msptsv_file_path + '.vit')
            vit_or_fbk_or_fbtsv_or_msptsv_list.append(1)
        elif path.exists(vit_fbk_fbtsv_msptsv_file_path + '.fbk'):
            vit_fbk_fbtsv_msptsv_file_list.append(vit_fbk_fbtsv_msptsv_file_path + '.fbk')
            vit_or_fbk_or_fbtsv_or_msptsv_list.append(2)
        elif path.exists(vit_fbk_fbtsv_msptsv_file_path + '.fb.tsv'):
            vit_fbk_fbtsv_msptsv_file_list.append(vit_fbk_fbtsv_msptsv_file_path + '.fb.tsv')
            vit_or_fbk_or_fbtsv_or_msptsv_list.append(3)
        elif path.exists(vit_fbk_fbtsv_msptsv_file_path + '.msp.tsv'):
            vit_fbk_fbtsv_msptsv_file_list.append(vit_fbk_fbtsv_msptsv_file_path + '.msp.tsv')
            vit_or_fbk_or_fbtsv_or_msptsv_list.append(4)
        else:
            vit_fbk_fbtsv_msptsv_file_list.append('')
            vit_or_fbk_or_fbtsv_or_msptsv_list.append(0)
    if (1 in beagle_or_vcf_list) and (2 in beagle_or_vcf_list):
        is_mixed = True
    else:
        is_mixed = False

    # Initialization:
    rs_ID_dict = {}
    masks =[]
    rs_ID_list = []
    ind_ID_list = []

    for i in range(num_arrays):
        logging.info("------ Array "+ str(i+1) + " Processing: ------")
        genome_matrix, ind_IDs, rs_IDs, rs_ID_dict = get_masked_matrix(beagle_vcf_file_list[i], beagle_or_vcf_list[i],
                                                                       vit_fbk_fbtsv_msptsv_file_list[i],
                                                                       vit_or_fbk_or_fbtsv_or_msptsv_list[i], is_mixed, is_masked,
                                                                       num_ancestries, average_parents, prob_thresh, rs_ID_dict,
                                                                       rsid_or_chrompos)


        masks.append(genome_matrix)
        rs_ID_list.append(rs_IDs)
        if (average_parents == False):
            ind_ID_list.append(ind_IDs)
        else:
            ind_ID_list.append(remove_AB_indIDs(ind_IDs))
        
    return masks, rs_ID_list, ind_ID_list

###########################################################################
def remove_AB_indIDs(ind_IDs):
    new_ind_IDs = []
    for i in range(int(len(ind_IDs)/2)):
        new_ind_IDs.append(ind_IDs[2*i][:-2])
    new_ind_IDs = np.array(new_ind_IDs)
    return new_ind_IDs

def add_AB_indIDs(ind_IDs):
    new_ind_IDs = []
    for i in range(len(ind_IDs)):
        new_ind_IDs.append(str(ind_IDs[i]) + '_A')
        new_ind_IDs.append(str(ind_IDs[i]) + '_B')
    new_ind_IDs = np.array(new_ind_IDs)
    return new_ind_IDs

def process_labels_weights(labels_file, masks, rs_ID_list, ind_ID_list, average_parents, num_arrays, ancestry, min_percent_snps, remove_labels_dict, is_weighted, save_masks, masks_file):
    labels_df = pd.read_csv(labels_file, sep='\t')
    label_list = []
    weight_list = []
    for array_ind in range(num_arrays):
        masked_matrix = masks[array_ind][ancestry]
        ind_IDs = ind_ID_list[array_ind]
        if average_parents:
            labels = np.array(labels_df['label'][labels_df['indID'].isin(ind_IDs)])
            label_ind_IDs = np.array(labels_df['indID'][labels_df['indID'].isin(ind_IDs)])
        else:
            temp_ind_IDs = remove_AB_indIDs(ind_IDs)
            labels = np.array(labels_df['label'][labels_df['indID'].isin(temp_ind_IDs)])
            labels = np.repeat(labels, 2)
            label_ind_IDs = np.array(labels_df['indID'][labels_df['indID'].isin(temp_ind_IDs)])
            label_ind_IDs = add_AB_indIDs(label_ind_IDs)
        keep_indices = [ind_IDs.tolist().index(x) for x in label_ind_IDs]
        masked_matrix = masked_matrix[:,keep_indices]
        ind_IDs = ind_IDs[keep_indices]
        array_num = array_ind + 1
        if not is_weighted:
            weights = np.ones(len(labels))
            combinations = np.zeros(len(labels))
            combination_weights = np.zeros(len(labels))
        else:
            if average_parents:
                weights = np.array(labels_df['weight'][labels_df['indID'].isin(ind_IDs)])
                if 'combination' in labels_df.columns:
                    combinations = np.array(labels_df['combination'][labels_df['indID'].isin(ind_IDs)])
                else:
                    combinations = np.zeros(len(weights))
                if 'combination_weight' in labels_df.columns:
                    combination_weights = np.array(labels_df['combination_weight'][labels_df['indID'].isin(ind_IDs)])
                else:
                    combination_weights = np.ones(len(weights))
            else:
                temp_ind_IDs = remove_AB_indIDs(ind_IDs)
                weights = np.array(labels_df['weight'][labels_df['indID'].isin(temp_ind_IDs)])
                weights = np.repeat(weights, 2)
                if 'combination' in labels_df.columns:
                    combinations = np.array(labels_df['combination'][labels_df['indID'].isin(temp_ind_IDs)])
                    combinations = np.repeat(combinations, 2)
                else:
                    combinations = np.zeros(len(weights))
                if 'combination_weight' in labels_df.columns:
                    combination_weights = np.array(labels_df['combination_weight'][labels_df['indID'].isin(temp_ind_IDs)])
                    combination_weights = np.repeat(combination_weights, 2)
                else:
                    combination_weights = np.ones(len(weights))
        if array_num in remove_labels_dict:
            remove_labels = remove_labels_dict[array_num]
            for i in range(len(labels)):
                if labels[i] in remove_labels:
                    weights[i] = 0
        percent_snps = 100 * (1 - np.mean(np.isnan(masked_matrix), axis=0))
        keep_indices = np.argwhere(percent_snps >= min_percent_snps).flatten()
        masked_matrix = masked_matrix[:,keep_indices]
        ind_IDs = ind_IDs[keep_indices]
        labels = labels[keep_indices]
        weights = weights[keep_indices]
        combinations = combinations[keep_indices]
        combination_weights = combination_weights[keep_indices]
        keep_indices = np.argwhere(weights > 0).flatten()
        masked_matrix_new = masked_matrix[:,keep_indices]
        ind_IDs_new = ind_IDs[keep_indices]
        labels_new = labels[keep_indices]
        weights_new = weights[keep_indices]
        pos_combinations = sorted(set(combinations[combinations > 0]))
        num_combinations = len(pos_combinations)
        if num_combinations > 0:
            for combination in pos_combinations:
                combined_indices = np.argwhere(combinations == combination)
                combined_col = np.nanmean(masked_matrix[:,combined_indices], axis=1)
                masked_matrix_new = np.append(masked_matrix_new, combined_col, axis=1)
                ind_IDs_new = np.append(ind_IDs_new, 'combined_ind_' + str(combination))
                labels_new = np.append(labels_new, labels[combined_indices[0][0]])
                weights_new = np.append(weights_new, combination_weights[combined_indices[0][0]])
        masked_matrix = masked_matrix_new
        ind_IDs = ind_IDs_new
        labels = labels_new
        weights = weights_new
        masks[array_ind][ancestry] = masked_matrix
        ind_ID_list[array_ind] = ind_IDs
        label_list += labels.tolist()
        weight_list += weights.tolist()
    label_list = np.array(label_list)
    weight_list = np.array(weight_list)
    if save_masks:
        np.savez_compressed(masks_file, masks=masks, rs_ID_list=rs_ID_list, ind_ID_list=ind_ID_list,
                 labels=label_list, weights=weight_list)
    return masks, ind_ID_list, label_list, weight_list

def center_masked_matrix(masked_matrix):
    masked_matrix -= np.nanmean(masked_matrix, axis=0)
    return masked_matrix

###########################################################################
def logger_config(verbose=True):
    logging_config = {"version": 1, "disable_existing_loggers": False}
    fmt = '[%(levelname)s] %(asctime)s: %(message)s'
    logging_config["formatters"] = {"basic": {"format": fmt, "datefmt": "%Y-%m-%d %H:%M:%S"}}
    now = datetime.datetime.now()

    logging_config["handlers"] = {
            "console": {
                "class": "logging.StreamHandler",
                "level": "DEBUG" if verbose else "INFO",
                "formatter": "basic",
                "stream": "ext://sys.stdout"
            },
            "info_file_handler": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "DEBUG" if verbose else "INFO", 
                "formatter": "basic",
                "maxBytes": 10485760,
                "backupCount": 20,
                "filename": f"log{now.year}_{now.month}_{now.day}__{now.hour}_{now.minute}.txt", # choose a better name or name as param?
                "encoding": "utf8"
                }
            }
    logging_config["root"] = {
                "level": "DEBUG",
                "handlers": ["console", "info_file_handler"]
            }
    return logging_config