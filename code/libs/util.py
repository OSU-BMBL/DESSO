from __future__ import division
import sys
import os
import gzip
import csv
import math
import numpy as np
from constants import *
import tensorflow as tf
import collections
from Bio import SeqIO
from Bio import SeqUtils
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC
from Bio.SubsMat import FreqTable

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


# Make dir if the given path is not exist
def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)


# Get ChIP-seq data list
def get_data_name():
    train_data_list = []
    table_name_list = []
    with open(PATH_ENCODE_TFBS) as f:
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

    for chr_index in chr_set:
        chrome[chr_index]['start'] = np.array(chrome[chr_index]['start'])
        chrome[chr_index]['end'] = np.array(chrome[chr_index]['end'])
        
    return chrome, peak_num


# Load sequences, targets, and DNA shapes from ENCODE
def load_data_encode(data_path, peak_coor, data_name, path_curr_data, train_test, dna_shape_flag, peak_flank, back_grou):
    # Make dir for training or test data
    make_dir(path_curr_data + "/" + train_test)
    peak_length = 2 * peak_flank + 1

    # Check whether the sequences.npy (contains combination of the training and validation data) already exists.
    # If no, read and process the peak sequences, otherwise, read sequences.npy, targets.npy, and sequences_alph.npy directly.
    if not os.path.exists(path_curr_data + "/" + train_test + "/" + "sequences.npy"):
        # Read the sequence file
        with gzip.open(data_path, 'r') as f:
            sequenceList = list(csv.reader(f, delimiter = '\t')) 
            invalid_seq_index = [i for i, seq in enumerate(sequenceList[1:]) if 'N' in seq[2].strip()]  

            # If the number of valid sequences less than 10000, we generate (10000 - #valid sequences) sequences by resampling from the available sequences
            if train_test == "train" and len(sequenceList) - 1 - len(invalid_seq_index) < MINIMUM_POSI_SEQ:
                sample_index = np.random.choice(len(sequenceList) - 1 - len(invalid_seq_index), MINIMUM_POSI_SEQ - (len(sequenceList) - 1 - len(invalid_seq_index)), replace = True)
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
            sequences = np.array([[float(x) for x in seq] for seq in sequences])    

            # Remove invalid targets according to invalid_seq_index
            targets = np.delete(targets, invalid_seq_index, 0) 

        # Generate negative sequences
        if train_test != 'all':
            seqs = []
            while len(seqs) < sequences.shape[0]:
                GC_content = SeqUtils.GC(sequences_alph[len(seqs)])                        # GC content of current DNA sequence 
                chr_index = chrom_seqSize.keys()[np.random.randint(len(chrom_seqSize))]    # Randomly select a chromosome
                chr_start = np.random.randint(chrom_seqSize[chr_index] - peak_length)      # Randomly select a start position in the selected chromosome
                chr_end = chr_start + peak_length                                          # The end position of the selected segment
                seq_segment = str(genome_seq[chr_index][chr_start : chr_end]).upper()

                # The sequences without overlapping and only consisting of 'ACGT' will be used as background sequence                    
                if ((chr_index not in peak_coor.keys()) or (not (np.any((peak_coor[chr_index]['start'] <= chr_start) & (chr_start <= peak_coor[chr_index]['end'])) or np.any((peak_coor[chr_index]['start'] <= chr_end) & (chr_end <= peak_coor[chr_index]['end'])) or np.any((peak_coor[chr_index]['start'] >= chr_start) & (peak_coor[chr_index]['end'] <= chr_end))))) and (len(set(seq_segment) - set('ACGT')) == 0) and abs(SeqUtils.GC(seq_segment) - GC_content) < 1:
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
        # Using Feature_encoding to generate DNA shape features    
        fdream5 = open(path_curr_data + "/" + train_test + '/encode.txt.s', 'w')
        for i in range(len(sequences_alph)):
            fdream5.write(str(1) + '  ' + sequences_alph[i] + '\n')
        fdream5.close()

        curr_path = os.getcwd()
        os.chdir(path_curr_data + "/" + train_test)
        os.system(PATH_LIBS + "/Feature_encoding/Feature_encoding.pl " + path_curr_data + "/" + train_test + "/encode.txt.s")
        os.chdir(curr_path)
        
    # Extract DNA shape features, including MGW, ProT, Roll, and HelT
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
    else:
        HelT = MGW = ProT = Roll = dnaShape = None
    return sequences, targets, sequences_alph, HelT, MGW, ProT, Roll, dnaShape


# Add appendix to the given data
def add_appendix(data, appendix_base, FILTER_LENGTH):
    new_data = np.ndarray(shape = (data.shape[0], data.shape[1] + 2 * (FILTER_LENGTH - 1), data.shape[2], data.shape[3]), dtype = np.float32)
    appendix = appendix_base * np.ones((FILTER_LENGTH - 1, data.shape[2]))

    for i in range(data.shape[3]):
        new_data[..., i] = np.array([np.concatenate((appendix, seq_array, appendix), axis = 0) for seq_array in data[..., i]])
    return new_data


# Prepare training data and test data based on feature format
def prep_data_encode(data_array, feature_format):
    data_comp = {}
    if feature_format == 'Seq':
        data_comp['0'] = data_array['Seq']
    elif feature_format == 'DNAShape':
        data_comp['0'] = data_array['HelT']
        data_comp['1'] = data_array['MGW']
        data_comp['2'] = data_array['ProT']
        data_comp['3'] = data_array['Roll']
    return data_comp


# Convert sequences to one-hot format
def oneHot_data_encode(sequences):
    return (np.arange(1, SEQUENCE_WIDTH + 1) == sequences.flatten()[:, None]).astype(np.float32).reshape(len(sequences), sequences.shape[1], SEQUENCE_WIDTH)[..., None]


def para_trainable(data_shape, init_sd, initValue, scope_name):
    with tf.variable_scope(scope_name, reuse=True):
        weights = tf.Variable(tf.truncated_normal(data_shape, stddev = init_sd, dtype = tf.float32, seed = SEED), name = 'weights')
        biases = tf.Variable(tf.constant(initValue, shape = [data_shape[-1]], dtype = tf.float32), name = 'biases')
    return weights, biases


# Convert normalized DNA shape features to their original scale
def convert_DNAShape(normalized_DNAShape, feature_format, key):
    if feature_format == "HelT" or (feature_format == "DNAShape" and key == "0"):
        return (normalized_DNAShape * HelT_span_first + HelT_min_first)
    elif feature_format == "MGW" or (feature_format == "DNAShape" and key == "1"):
        return (normalized_DNAShape * MGW_span_first + MGW_min_first)
    elif feature_format == "ProT" or (feature_format == "DNAShape" and key == "2"):
        return (normalized_DNAShape * ProT_span_first + ProT_min_first)
    elif feature_format == "Roll" or (feature_format == "DNAShape" and key == "3"):
        return (normalized_DNAShape * Roll_span_first + Roll_min_first)


# Calculate infomration content of given sequences
def cal_IC(path_motifInstance):    
    e_freq_table = FreqTable.FreqTable(EXPECT_FREQ, FreqTable.FREQ, IUPAC.unambiguous_dna)
    information_content = []
    alignment = AlignIO.read(path_motifInstance, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)

    for j in range(FILTER_LENGTH):                
        information_content.append(summary_align.information_content(j, j + 1, e_freq_table = e_freq_table, chars_to_ignore = ['N']))
    return information_content
