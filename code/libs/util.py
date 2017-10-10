### Some utils used in other scripts
from __future__ import division
import sys
import re
import os
import gzip
import csv
import time
import math
import numpy as np
import subprocess
import cPickle
from global_cons import *
from dream5_cons import *
from DNAShape_cons import *
#from serie2image import *
import encode_cons as encode
import sequence
import tensorflow as tf
import collections
sys.path.append(PATH_LIBS + '/biopython')
from Bio import SeqIO

########################################################
####################### ENCODE #########################
########################################################
# Get human genome sequence and shape
path_genome = PATH_DATA + '/GRCh37.p13.genome.fa'    # Fasta files of Release 19 (GRCh37.p13)
genome_seq = {}                                      # Human genome sequence, where key = chromosome id, value = sequence
genome_shape = collections.defaultdict(dict)         # Human genome shape, where key = chromosome id, value = DNA shape
chrom_seqSize = {}                                   # Size of each chromosome, where key = chromosome id, value = length

# Read genome sequence from human genome
fasta_sequences = SeqIO.parse(open(path_genome),'fasta')
for fasta in fasta_sequences:
    if fasta.id.startswith('chr'):
        chrom_seqSize[fasta.id] = len(fasta)
        genome_seq[fasta.id] = fasta.seq
'''
start_time = time.time()
# Read DNA shape of human genome
for chr_index in chrom_seqSize.keys():
    for DNAshape in encode.DNASHAPE_FEATURE:
        genome_shape[chr_index][DNAshape] = np.load(PATH_DATA + "/GRCh37.p13.genome/" + chr_index + "/" + chr_index + "." + DNAshape + ".npy")
print('Time for DNA shape load: %d seconds' % (time.time() - start_time))
'''
# Make dir if the given path is not exist
def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)


# Get ChIP-seq data list
def get_data_name():
    train_data_list = []
    table_name_list = []
    with open(encode.PATH_ENCODE_TFBS) as f:
        results = list(csv.reader(f, delimiter = '\t'))
        
    for i in range(len(results)):
        train_data_list.append(results[i][0].strip())
        table_name_list.append(results[i][1].strip())
    return train_data_list, table_name_list


# Coordinates of each peak
def get_peak_coor(tfbs_unif_path):
    # Get all chromosome included in this file 
    with gzip.open(tfbs_unif_path, 'r') as fnarrow_peak:
        peak_info = list(csv.reader(fnarrow_peak, delimiter = '\t'))
    peak_num = len(peak_info)    # Number of peaks
    chr_set = set([peak_info[i][0] for i in range(len(peak_info))])
            
    # Construct a nested dictionalry to record gene position
    chrome = collections.defaultdict(dict)
    for chr_index in chr_set:
        chrome[chr_index]['start'] = []
        chrome[chr_index]['end'] = []

    # Record coordinates of each peak  
    for line in peak_info:
        chr_start = line[1]         # Start position of current peak
        chr_end = line[2]           # End position of current peak
    
        chrome[line[0]]['start'].append(int(chr_start))
        chrome[line[0]]['end'].append(int(chr_end))
    return chrome, peak_num


# Load sequences, targets, and DNA shapes from ENCODE
def load_data_encode(data_path, peak_coor, data_name, path_curr_data, train_test, dna_shape_flag, peak_flank, back_grou):
    # Make dir for training or testing data
    make_dir(path_curr_data + "/" + train_test)
    peak_length = 2 * peak_flank + 1

    # Check whether the sequences.npy (contains combination of the training and validation data) already exists.
    # If no, read and process the peak sequences, otherwise, read sequences.npy, targets.npy, and sequences_alph.npy directly.
    # These three npy files are synchronous, so here we only check the existence of the sequences.npy
    if not os.path.exists(path_curr_data + "/" + train_test + "/" + "sequences.npy"):
        # Read the sequence file
        with gzip.open(data_path, 'r') as f:
            sequenceList = list(csv.reader(f, delimiter = '\t')) 
            invalid_seq_index = [i for i, seq in enumerate(sequenceList[1:]) if 'N' in seq[2].strip()]    # Some sequences may have 'N'

            # If the number of valid sequences less than 10000, we generate (10000 - #valid sequences) sequences by resampling from the available sequences
            if train_test == "train" and len(sequenceList) - 1 - len(invalid_seq_index) < encode.MINIMUM_POSI_SEQ:
                sample_index = np.random.choice(len(sequenceList) - 1 - len(invalid_seq_index), encode.MINIMUM_POSI_SEQ - (len(sequenceList) - 1 - len(invalid_seq_index)), replace = True)
                sequenceList_valid = [sequenceList[1:][i] for i in xrange(len(sequenceList[1:])) if i not in invalid_seq_index]
                sample_sequence = [sequenceList_valid[sample_index[i]] for i in xrange(len(sample_index))]
                sequenceList = sequenceList + sample_sequence
            
            # Sequences (numpy array)
            sequences = np.array([[float(BASE_NUM[x]) for x in sequenceList[i][2]] for i in xrange(1, len(sequenceList))])

            # Sequences ('A', 'T', 'C', 'G')
            sequences_alph = [sequenceList[i][2] for i in xrange(1, len(sequenceList))]

            # Targets (0/1)
            targets = np.array([int(sequenceList[i][3]) for i in xrange(1, len(sequenceList))])
            
            # Remove invalid sequences according to invalid_seq_index
            sequences_alph = np.array(np.delete(sequences_alph, invalid_seq_index, 0))
                       
            # Remove invalid sequences according to invalid_seq_index
            sequences = np.array(np.delete(sequences, invalid_seq_index, 0))
            sequences = np.array([[float(x) for x in seq] for seq in sequences])    # [type: numpy.ndarray]

            # Remove invalid targets according to invalid_seq_index
            targets = np.delete(targets, invalid_seq_index, 0) 

        if back_grou == 'dinu_shuf' and train_test != 'all':
            #####################################################################################
            ### 1. Generate negative sequences using dinucleotide-preserving shuffle strategy ###        
            f_fa = open(path_curr_data + "/" + train_test + '/seq.fa', 'w')
            for i in range(len(sequences_alph)):
                f_fa.write('>Seq' + str(i) + '\n')
                f_fa.write(sequences_alph[i] + '\n')
            f_fa.close()

            # Shuffle the sequences using fasta-dinucleotide-shuffle.py which is a tool downloaded from MEME
            featEnco = subprocess.Popen([PATH_LIBS + "/fasta-dinucleotide-shuffle.py", "-f",  path_curr_data + "/" + train_test + "/seq.fa", "-o", path_curr_data + "/" + train_test + "/nega_seq.fa"])
            featEnco.wait()

            # Then conbine the original data (positive class) and the shuffled data (negative class) together.
            # Firstly, read sequences using "sequence" which is a module downloaded from MEME
            seqs = sequence.readFASTA(path_curr_data + "/" + train_test + '/nega_seq.fa', 'Extended DNA')
            
            seq_shuffle = np.array([[float(BASE_NUM[x]) for x in seq.getString()] for seq in seqs])
            seq_alph_shuffle = [seq.getString() for seq in seqs]
            
            # Concatenate the original sequences and the shuffled sequences
            sequences = np.concatenate((sequences, seq_shuffle), axis = 0)

            # Concatenate the original targets and the all-zero targets
            targets = np.concatenate((targets, np.zeros((seq_shuffle.shape[0],), dtype = np.int)), axis = 0)

            # Concatenate the original sequences_alph and the shuffled sequences_alph
            sequences_alph = np.concatenate((sequences_alph, seq_alph_shuffle), axis = 0)
        
        elif back_grou == 'rand_geno' and train_test != 'all':
            seqs = []    # Store the sequences which are generated by randomly selecting (2 * PEAK_FLANK + 1) bps in human genome
            # These generate sequences will be used as background sequences in binary classification;
            # Here, we generate the same number of background sequences as foreground sequences
            while len(seqs) < sequences.shape[0]:
                chr_index = chrom_seqSize.keys()[np.random.randint(len(chrom_seqSize))]                  # Randomly select a chromosome
                chr_start = np.random.randint(chrom_seqSize[chr_index] - peak_length)    # Randomly select a start position in selected chromosome
                chr_end = chr_start + peak_length                                        # The end position of the selected segment
                seq_segment = str(genome_seq[chr_index][chr_start : chr_end]).upper()

                # The sequences without overlapping and only consisting of 'ACGT' will be used as background sequence                    
                if ((chr_index not in peak_coor.keys()) or (not (np.any((peak_coor[chr_index]['start'] <= chr_start) & (chr_start <= peak_coor[chr_index]['end'])) or np.any((peak_coor[chr_index]['start'] <= chr_end) & (chr_end <= peak_coor[chr_index]['end'])) or np.any((peak_coor[chr_index]['start'] >= chr_start) & (peak_coor[chr_index]['end'] <= chr_end))))) and (len(set(seq_segment) - set('ACGT')) == 0):
                    seqs.append(seq_segment)
                            
            seq_shuffle = np.array([[float(BASE_NUM[x]) for x in seq] for seq in seqs])
            seq_alph_shuffle = [seq for seq in seqs]

            # Concatenate the original sequences and the shuffled sequences
            sequences = np.concatenate((sequences, seq_shuffle), axis = 0)

            # Concatenate the original targets and the all-zero targets
            targets = np.concatenate((targets, np.zeros((seq_shuffle.shape[0],), dtype = np.int)), axis = 0)

            # Concatenate the original sequences_alph and the shuffled sequences_alph
            sequences_alph = np.concatenate((sequences_alph, seq_alph_shuffle), axis = 0)

        np.save(path_curr_data + "/" + train_test + "/" + "sequences.npy", sequences)
        np.save(path_curr_data + "/" + train_test + "/" + "targets.npy", targets)
        np.save(path_curr_data + "/" + train_test + "/" + "sequences_alph.npy", sequences_alph)

    # Obtain the data if they have already been saved
    sequences = np.load(path_curr_data + "/" + train_test + "/" + "sequences.npy")
    targets = np.load(path_curr_data + "/" + train_test + "/" + "targets.npy")
    sequences_alph = np.load(path_curr_data + "/" + train_test + "/" + "sequences_alph.npy")
    
    if dna_shape_flag and not os.path.exists(path_curr_data + "/" + train_test + "/encode_4shape.matlab"):
        ########### DNA Shape from Feature_encoding ##############
        ##########################################################
        # Obtain the DNA shape
        # Using Feature_encoding to generate the DNA shape       
        fdream5 = open(path_curr_data + "/" + train_test + '/encode.txt.s', 'w')
        for i in range(len(sequences_alph)):
            fdream5.write(str(1) + '  ' + sequences_alph[i] + '\n')
        fdream5.close()
        featEnco = subprocess.Popen([PATH_LIBS + "/Feature_encoding/Feature_encoding.pl", path_curr_data + "/" + train_test + "/encode.txt.s", ">", path_curr_data + "/" + train_test])
        featEnco.wait()
        
    # Extract the DNA shape, including MGW, ProT, Roll, and HelT
    if dna_shape_flag:
        with open(path_curr_data + "/" + train_test + '/encode_4shape.matlab') as f:
            reader = csv.reader(f, delimiter = ' ')
            dnaShape = np.array([[float(x) for x in row] for row in reader])  
    
        # The first column are all-one, so delete them firstly
        dnaShape = np.delete(dnaShape, 0, 1)

        MGW = dnaShape[:, 0 : peak_length - 4]
        Roll = dnaShape[:, peak_length - 4 : (2 * peak_length) - 7]
        HelT = dnaShape[:, (2 * peak_length) - 7 : (3 * peak_length) - 10]
        ProT = dnaShape[:, (3 * peak_length) - 10 : (4 * peak_length) - 14]

        MGW = MGW.reshape(MGW.shape[0], MGW.shape[1], 1)
        Roll = Roll.reshape(Roll.shape[0], Roll.shape[1], 1)
        HelT = HelT.reshape(HelT.shape[0], HelT.shape[1], 1)
        ProT = ProT.reshape(ProT.shape[0], ProT.shape[1], 1)

        #MGW_2nd = np.multiply(MGW[:, 0 : MGW.shape[1] - 1], MGW[:, 1 : MGW.shape[1]])
        #Roll_2nd = np.multiply(Roll[:, 0 : Roll.shape[1] - 1], Roll[:, 1 : Roll.shape[1]])
        #HelT_2nd = np.multiply(HelT[:, 0 : HelT.shape[1] - 1], HelT[:, 1 : HelT.shape[1]])
        #ProT_2nd = np.multiply(ProT[:, 0 : ProT.shape[1] - 1], ProT[:, 1 : ProT.shape[1]])
        
        #sequences_array = (np.arange(1, SEQUENCE_WIDTH + 1) == sequences.flatten()[:, None]).astype(np.float32).reshape(len(sequences), sequences.shape[1] * SEQUENCE_WIDTH)
        #dnaShape = dnaShape_image
        #dnaShape = sequences_array
        #dnaShape = np.concatenate((MGW_2nd, Roll_2nd, HelT_2nd, ProT_2nd), axis = 1)
        #dnaShape = np.concatenate((dnaShape, MGW_2nd, Roll_2nd, HelT_2nd, ProT_2nd), axis = 1)
        #dnaShape = np.concatenate((sequences_array, dnaShape, MGW_2nd, Roll_2nd, HelT_2nd, ProT_2nd), axis = 1)
    else:
        HelT = MGW = ProT = Roll = dnaShape = None

    return sequences, targets, sequences_alph, HelT, MGW, ProT, Roll, dnaShape


# Add appendix to the data
def add_appendix(data, appendix_base, filter_length):
    print(data.shape)
    new_data = np.ndarray(shape = (data.shape[0], data.shape[1] + 2 * (filter_length - 1), data.shape[2], data.shape[3]), dtype = np.float32)
    appendix = appendix_base * np.ones((filter_length - 1, data.shape[2]))

    for i in range(data.shape[3]):
        new_data[..., i] = np.array([np.concatenate((appendix, seq_array, appendix), axis = 0) for seq_array in data[..., i]])
    return new_data


# Extract training data, training targets, validation data, validation targets, test data, and test targets
def extract_data_encode(feature_format, train_sequences_array, train_MGW_array, train_ProT_array, train_targets_array, train_dnaShape, train_MP_array, test_sequences_array, test_MGW_array, test_ProT_array, test_targets_array, test_dnaShape, test_MP_array, dna_shape_flag):
    random_index = np.random.permutation(train_sequences_array.shape[0])
    
    sequences_array_train = train_sequences_array[random_index, ...]
    sequences_array_test = test_sequences_array

    targets_array_train = train_targets_array[random_index, ...]
    targets_array_test = test_targets_array

    if dna_shape_flag:
        MGW_array_train = train_MGW_array[random_index, ...]
        MGW_array_test = test_MGW_array
            
        ProT_array_train = train_ProT_array[random_index, ...]
        ProT_array_test = test_ProT_array
            
        dnaShape_array_train = train_dnaShape[random_index, ...]
        dnaShape_array_test = test_dnaShape

        MP_array_train = train_MP_array[random_index, ...]
        MP_array_test = test_MP_array
                
    # The feature_format should be like 'S' (sequences), 'M' (MGW), 'P' (ProT), 'SM' (sequences and MGW), 'SP' (sequences and ProT), 'MP' (MGW and ProT), 'SMP' (sequences, MGW and ProT)
    if (len(feature_format) == 1):
        if (feature_format == 'S'):    # if we only use sequences as features
            
            train_data = sequences_array_train[encode.VALIDATION_SIZE:, ..., None]
            train_targets =  targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = sequences_array_train[:encode.VALIDATION_SIZE, ..., None]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = sequences_array_test[..., None]
            test_targets = targets_array_test

            if dna_shape_flag:
                train_dnaShape = dnaShape_array_train[encode.VALIDATION_SIZE:, ..., None]
                validation_dnaShape = dnaShape_array_train[:encode.VALIDATION_SIZE, ..., None]
                test_dnaShape = dnaShape_array_test[..., None]
                '''
                [S, V] = PCA(train_dnaShape)
                K = 500
                train_dnaShape = proj_data(train_dnaShape, V, K)
                validation_dnaShape = proj_data(validation_dnaShape, V, K)
                test_dnaShape = proj_data(test_dnaShape, V, K)
                '''
                train_MP = MP_array_train[encode.VALIDATION_SIZE:, ..., None]
                validation_MP = MP_array_train[:encode.VALIDATION_SIZE, ..., None]
                test_MP = MP_array_test[..., None]
            else:
                train_dnaShape = validation_dnaShape = test_dnaShape = train_MP = validation_MP = test_MP = None
            '''
            train_data = sequences_array_train[VALIDATION_SIZE:, 2 : VALID_SEQUENCE_LENGTH - 2, ..., None]
            train_targets =  targets_array_train[VALIDATION_SIZE:, ...]
            validation_data = sequences_array_train[:VALIDATION_SIZE, 2 : VALID_SEQUENCE_LENGTH - 2, ..., None]
            validation_targets = targets_array_train[:VALIDATION_SIZE, ...]
            test_data = sequences_array_test[:, 2 : VALID_SEQUENCE_LENGTH - 2, ..., None]
            test_targets = targets_array_test
            '''
        elif (feature_format == 'M'):    # if we only use MGW as features
            train_data = MGW_array_train[encode.VALIDATION_SIZE:, ..., None]
            train_targets = targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = MGW_array_train[:encode.VALIDATION_SIZE, ..., None]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = MGW_array_test[..., None]
            test_targets = targets_array_test
            
        elif (feature_format == 'P'):    # if we only use ProT as features
            train_data = ProT_array_train[encode.VALIDATION_SIZE:, ..., None]
            train_targets = targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = ProT_array_train[:encode.VALIDATION_SIZE, ..., None]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = ProT_array_test[..., None]
            test_targets = targets_array_test    
    else:
        # Stack sequences, MGW, and ProT together
        train_validation_data = np.ndarray(shape = np.append(MGW_array_train.shape, 3), dtype = np.float32)
        train_validation_data[..., 0] = sequences_array_train[:, 2:encode.VALID_SEQUENCE_LENGTH - 2, ...]
        train_validation_data[..., 1] = MGW_array_train
        train_validation_data[..., 2] = ProT_array_train
        
        test_stack = np.ndarray(shape = np.append(MGW_array_test.shape, 3), dtype = np.float32)
        test_stack[..., 0] = sequences_array_test[:, 2:encode.VALID_SEQUENCE_LENGTH - 2, ...]
        test_stack[..., 1] = MGW_array_test
        test_stack[..., 2] = ProT_array_test
        
        if (feature_format == 'SM'):    # If we use sequences and MGW as features
            train_data = train_validation_data[encode.VALIDATION_SIZE:, ..., [0, 1]]
            train_targets = targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = train_validation_data[:encode.VALIDATION_SIZE, ..., [0, 1]]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = test_stack[:, ..., [0, 1]]
            test_targets = targets_array_test
            
        elif (feature_format == 'SP'):    # If we use sequences and ProT as features
            train_data = train_validation_data[encode.VALIDATION_SIZE:, ..., [0, 2]]
            train_targets = targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = train_validation_data[:encode.VALIDATION_SIZE, ..., [0, 2]]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = test_stack[:, ..., [0, 2]]
            test_targets = targets_array_test
            
        elif (feature_format == 'MP'):    # If we use ProT and MGW as features
            train_data = train_validation_data[encode.VALIDATION_SIZE:, ..., [1, 2]]
            train_targets = targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = train_validation_data[:encode.VALIDATION_SIZE, ..., [1, 2]]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = test_stack[:, ..., [1, 2]]
            test_targets = targets_array_test
            
        elif (feature_format == 'SMP'):    # If we use sequences, MGW and ProT as features
            train_data = train_validation_data[encode.VALIDATION_SIZE:, ..., [0, 1, 2]]
            train_targets = targets_array_train[encode.VALIDATION_SIZE:, ...]
            validation_data = train_validation_data[:encode.VALIDATION_SIZE, ..., [0, 1, 2]]
            validation_targets = targets_array_train[:encode.VALIDATION_SIZE, ...]
            test_data = test_stack[:, ..., [0, 1, 2]]
            test_targets = targets_array_test
        
    ### Normalize the training data and test data
    #train_data, validation_data, test_data = normalize_x(train_data, validation_data, test_data)
    
    return train_data, train_targets, validation_data, validation_targets, test_data, test_targets, train_dnaShape, validation_dnaShape, test_dnaShape, train_MP, validation_MP, test_MP


# Prepare training data and test data based on feature format
# The feature format could be one of the ['Seq', 'HelT', 'MGW', 'ProT', 'Roll', 'DNAShape', 'Seq_DNAShape']
def prep_data_encode(data_array, feature_format):
    data_comp = {}
    if feature_format == 'Seq':
        data_comp['0'] = data_array['Seq']
    elif feature_format == 'HelT':
        data_comp['0'] = data_array['HelT']
    elif feature_format == 'MGW':
        data_comp['0'] = data_array['MGW']
    elif feature_format == 'ProT':
        data_comp['0'] = data_array['ProT']
    elif feature_format == 'Roll':
        data_comp['0'] = data_array['Roll']
    elif feature_format == 'DNAShape':
        data_comp['0'] = data_array['HelT']
        data_comp['1'] = data_array['MGW']
        data_comp['2'] = data_array['ProT']
        data_comp['3'] = data_array['Roll']
    elif feature_format == 'Seq_DNAShape':
        data_comp['0'] = data_array['Seq']
        data_comp['1'] = data_array['HelT']
        data_comp['2'] = data_array['MGW']
        data_comp['3'] = data_array['ProT']
        data_comp['4'] = data_array['Roll']
    return data_comp


# Convert the sequences and DNA shape (MGW, ProT) to one-hot format
def oneHot_data_encode(sequences, targets, HelT, MGW, ProT, Roll, dna_shape_flag):
    sequences_array = (np.arange(1, SEQUENCE_WIDTH + 1) == sequences.flatten()[:, None]).astype(np.float32).reshape(len(sequences), sequences.shape[1], SEQUENCE_WIDTH)
    sequences_array = sequences_array[..., None]
    targets_array = targets
    if dna_shape_flag:
        HelT_array = HelT[..., None]
        MGW_array = MGW[..., None]
        ProT_array = ProT[..., None]
        Roll_array = Roll[..., None]
        '''
        MGW_bina = (MGW > 0.5) * 1
        ProT_bina = (ProT > 0.5) * 1
        HelT_bina = (HelT > 0.5) * 1
        Roll_bina = (Roll > 0.5) * 1
        
        HelT_bina = np.concatenate((np.zeros(shape = (HelT_bina.shape[0], 1), dtype = np.int) ,HelT_bina), axis = 1)
        Roll_bina = np.concatenate((np.zeros(shape = (Roll_bina.shape[0], 1), dtype = np.int) ,Roll_bina), axis = 1)
        #sequences_array = np.array([np.concatenate((np.transpose(MGW_bina[i][None, :]), sequences_array[i], np.transpose(ProT_bina[i][None, :]), np.transpose(HelT_bina[i][None, :]), np.transpose(Roll_bina[i][None, :])), axis = 1) for i in xrange(len(sequences_array))])
        #sequences_array = np.array([np.concatenate((np.transpose(np.concatenate(([0], HelT[i]), axis = 0)[None, :]), sequences_array[i], np.transpose(np.concatenate(([0], Roll[i]), axis = 0)[None, :])), axis = 1) for i in xrange(len(sequences_array))])
        #sequences_array = np.array([np.concatenate((np.transpose(MGW[i][None, :]), np.transpose(np.concatenate(([0], HelT[i]), axis = 0)[None, :]), sequences_array[i], np.transpose(np.concatenate(([0], Roll[i]), axis = 0)[None, :]), np.transpose(ProT[i][None, :])), axis = 1) for i in xrange(len(sequences_array))])
        MGW_array = sequences_array * MGW[:, :, None]      # [len(sequences), VALID_SEQUENCE_LENGTH - 4, SEQUENCE_WIDTH]
        ProT_array = sequences_array * ProT[:, :, None]    # [len(sequences), VALID_SEQUENCE_LENGTH - 4, SEQUENCE_WIDTH]
        #MP = np.array([np.transpose(np.concatenate((MGW_bina[i,...][None, ...], ProT_bina[i, ...][None, ...]), axis = 0)) for i in xrange(len(MGW))])
        MP = np.array([np.transpose(np.concatenate((MGW_bina[i,...][None, ...], ProT_bina[i, ...][None, ...], HelT_bina[i, ...][None, ...], Roll_bina[i, ...][None, ...]), axis = 0)) for i in xrange(len(MGW))])
        '''
    else:
        HelT_array = MGW_array = ProT_array = Roll_array = None
        
    return sequences_array, targets_array, HelT_array, MGW_array, ProT_array, Roll_array

# Calculate error rate of binary classification
def error_rate(predictions, labels):
    return (100.0 * np.sum(np.argmax(predictions, 1) == labels) / predictions.shape[0])


def para_trainable(data_shape, init_sd, initValue, scope_name):
    # The variance of the neurals in the network is 0.1
    with tf.variable_scope(scope_name, reuse=True):
        weights = tf.Variable(tf.truncated_normal(data_shape, stddev = init_sd, dtype = tf.float32, seed = encode.SEED), name = 'weights')
        biases = tf.Variable(tf.constant(initValue, shape = [data_shape[-1]], dtype = tf.float32), name = 'biases')
    return weights, biases

# Initialize parameters in fully-connected neural networks
def para_trainable_nn(data_shape, init_sd, initValue, scope_name):
    with tf.variable_scope(scope_name, reuse=True):
        weights = tf.Variable(tf.random_uniform(data_shape, minval = -math.sqrt(6 / (data_shape[0] + data_shape[1])), maxval = math.sqrt(6 / (data_shape[0] + data_shape[1])), dtype = tf.float32, seed = encode.SEED), name = 'weights')
        #weights = tf.Variable(tf.random_uniform(data_shape, minval = -0.15, maxval = 0.15, dtype = tf.float32, seed = encode.SEED), name = 'weights')
        biases = tf.Variable(tf.constant(initValue, shape = [data_shape[-1]], dtype = tf.float32), name = 'biases')
    return weights, biases


### ROC-AUC calculation from DeepBind ###
def calc_auc(z, y, want_curve = False):
    assert len(z) == len(y)
    order = np.argsort(z, axis=0, kind="mergesort")[::-1].ravel()
    z = z[order]
    y = y[order]

    # Accumulate the true positives with decreasing threshold
    tpr = y.cumsum()
    fpr = 1 + np.arange(len(y)).astype(y.dtype) - tpr

    # If curve doesn't have point at x=0, explicitly insert one
    if fpr[0] != 0:
        tpr = np.r_[0,tpr]
        fpr = np.r_[0,fpr]

    # If one of the classes was empty, return NaN
    if fpr[-1] == 0 or tpr[-1] == 0:
        return (np.nan,None) if want_curve else np.nan

    # Convert sums to rates
    tpr = tpr / tpr[-1]
    fpr = fpr / fpr[-1]

    # Calculate area under the curve using trapezoidal rule
    auc = np.trapz(tpr,fpr,axis=0)

    # Done!
    if want_curve:
        curve = np.hstack([fpr,tpr])
        return auc,curve
    return auc

def placeholder_dnaShape_encode(dnaShape_length):
    train_data_node = tf.placeholder(tf.float32, shape = (encode.BATCH_SIZE, 64, 64, 1))
    eval_data_node = tf.placeholder(tf.float32, shape = (encode.EVAL_BATCH_SIZE, 64, 64, 1)) 
    return train_data_node, eval_data_node


def placeholder_MP_encode(sequence_length, num_channels):
    train_MP_node = tf.placeholder(tf.float32, shape = (encode.BATCH_SIZE, sequence_length, 4, num_channels))
    eval_MP_node = tf.placeholder(tf.float32, shape = (encode.EVAL_BATCH_SIZE, sequence_length, 4, num_channels))
    return train_MP_node, eval_MP_node
'''
def conv_op(fan_in, shape):
    W = tf.Variable(tf.truncated_normal(shape, stddev = 0.1, dtype = tf.float32, seed = encode.SEED))
    b = tf.Variable(tf.constant(0, shape = [shape[-1]], dtype = tf.float32))
    return tf.add(tf.nn.conv2d(fan_in, W, strides=[1,1,1,1], padding='SAME'), b)
'''


def conv_op(fan_in, shape, name):
    W = tf.get_variable("%s_W"%name, shape, tf.float32, tf.random_normal_initializer(0.0, 0.1))
    b = tf.get_variable("%s_b"%name, shape[-1], tf.float32, tf.constant_initializer(0))
    return tf.add(tf.nn.conv2d(fan_in, W, strides=[1,1,1,1], padding='SAME'), b)


# Convert normalized DNA shape features to their original scale
def convert_DNAShape(normalized_DNAShape, feature_format, key):
    if feature_format == "HelT" or (feature_format == "DNAShape" and key == "0") or (feature_format == "Seq_DNAShape" and key == "1"):
        return (normalized_DNAShape * HelT_span_first + HelT_min_first)
    elif feature_format == "MGW" or (feature_format == "DNAShape" and key == "1") or (feature_format == "Seq_DNAShape" and key == "2"):
        return (normalized_DNAShape * MGW_span_first + MGW_min_first)
    elif feature_format == "ProT" or (feature_format == "DNAShape" and key == "2") or (feature_format == "Seq_DNAShape" and key == "3"):
        return (normalized_DNAShape * ProT_span_first + ProT_min_first)
    elif feature_format == "Roll" or (feature_format == "DNAShape" and key == "3") or (feature_format == "Seq_DNAShape" and key == "4"):
        return (normalized_DNAShape * Roll_span_first + Roll_min_first)

    
