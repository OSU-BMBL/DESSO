from __future__ import division
import sys
import os
import gzip
import csv

from libs.constants import *
import libs.util as util

import argparse
import numpy as np
from six.moves import xrange
import tensorflow as tf
from scipy.stats import binom_test

parser = argparse.ArgumentParser()
parser.add_argument("--start_index", help = "Start index of the 690 ENCODE ChIP-seq datasets.", type = int)
parser.add_argument("--end_index", help = "END index of the 690 ENCODE ChIP-seq datasets.", type = int)
parser.add_argument("--peak_flank", help = "Number of flanking base pairs at each side of peak summit.", type = int)
parser.add_argument("--network", help = "Neural network used in model training.")
parser.add_argument("--feature_format", help = "Feature format of the input.")
parser.add_argument("--start_cutoff", help = "Start of the motif cutoff interval.")
parser.add_argument("--end_cutoff", help = "End of the motif cutoff interval.")
parser.add_argument("--step_cutoff", help = "Step of the motif cutoff interval.")
args = parser.parse_args()

PEAK_FLANK = [args.peak_flank]
NETWORK = [args.network]
FEATURE_FORMAT = [args.feature_format]

# Flag for generating DNA shape features
if args.feature_format == 'Seq':
    DNA_SHAPE = False
elif args.feature_format == 'DNAShape':
    DNA_SHAPE = True

# To-be-optimized motif cutoff
motif_cutoff = np.arange(float(args.start_cutoff), float(args.end_cutoff), float(args.step_cutoff))

for peak_flank in PEAK_FLANK:
    PATH_ENCODE = os.path.join(PATH_DATA, 'encode_' + str(2 * peak_flank + 1)) 
    for back_grou in BACK_GROU:
        PATH_OUTPUT_ENCODE = os.path.join(PATH_OUTPUT, 'encode_' + str(2 * peak_flank + 1), back_grou)
  
        [train_data_list, table_name_list] = util.get_data_name()    # The name list of the training datasets

        for train_data_name in train_data_list[args.start_index : args.end_index]:                               
            table_name = table_name_list[train_data_list.index(train_data_name)]
            [peak_coor, peak_num] = util.get_peak_coor(PATH_ENCODE_TFBS_UNIF + "/" + table_name + ".narrowPeak.gz")
                            
            path_curr_data = PATH_OUTPUT_ENCODE + '/' + train_data_name[-16:]
            path_background_seq = PATH_DATA + "/encode_" + str(2 * peak_flank + 1) + "_background"
            
            ################# Load Data ##################          
            [all_sequences, all_targets, all_sequences_alph, all_HelT, all_MGW, all_ProT, all_Roll, all_dnaShape] = \
            util.load_data_encode(PATH_ENCODE + "/" + train_data_name + ".seq.gz", peak_coor, train_data_name, path_curr_data, "all", DNA_SHAPE, peak_flank, back_grou)
            
            ####### Convert data to one-hot format #######
            all_sequences_array = util.oneHot_data_encode(all_sequences)
            all_data_array = {}
            all_data_array['Seq'] = all_sequences_array
            if DNA_SHAPE:
                all_data_array['HelT'] = all_HelT[..., None]
                all_data_array['MGW'] = all_MGW[..., None]
                all_data_array['ProT'] = all_ProT[..., None]
                all_data_array['Roll'] = all_Roll[..., None]
            
            motif_sequences_alph = all_sequences_alph[: (500 if all_sequences_array.shape[0] > 500 else all_sequences_array.shape[0]), ...]
            motif_data_array = {key : value[: (500 if value.shape[0] > 500 else value.shape[0]), ...] for key, value in all_data_array.items()}
            
            # Load background datasets including sequence and DNA shape
            back_sequences_alph = np.load(path_background_seq + "/seq_alph.npy")
            back_sequences_array = util.oneHot_data_encode(np.load(path_background_seq + "/seq.npy"))  
            background_data_array = {}
            background_data_array['Seq'] = back_sequences_array
            if DNA_SHAPE:
                background_data_array['HelT'] = np.load(path_background_seq + "/HelT.npy")[..., None][..., None]
                background_data_array['MGW'] = np.load(path_background_seq + "/MGW.npy")[..., None][..., None]
                background_data_array['ProT'] = np.load(path_background_seq + "/ProT.npy")[..., None][..., None]
                background_data_array['Roll'] = np.load(path_background_seq + "/Roll.npy")[..., None][..., None]                                    
                
            for feature_format in FEATURE_FORMAT:
                motif_data_comp = util.prep_data_encode(motif_data_array, feature_format)
                data_keys = motif_data_comp.keys()
                motif_data = tuple(motif_data_comp[key] for key in data_keys)
                
                background_data_comp = util.prep_data_encode(background_data_array, feature_format)
                background_data = tuple(background_data_comp[key] for key in data_keys)
                
                # Different neural network architectures
                for net in NETWORK:
                    path_curr_net = path_curr_data + "/" + feature_format + "/" + net                                            
                    file_name = "model.ckpt"
                    
                    with tf.Graph().as_default():
                        ### Placeholder inputs which will be fed into the graph ###
                        motif_data_node = tuple(tf.placeholder(tf.float32, shape = (motif_data_comp[key].shape[0], motif_data_comp[key].shape[1], motif_data_comp[key].shape[2], motif_data_comp[key].shape[3])) for key in data_keys)
                        motif_data_background_node = tuple(tf.placeholder(tf.float32, shape = (background_data_comp[key].shape[0], background_data_comp[key].shape[1], background_data_comp[key].shape[2], background_data_comp[key].shape[3])) for key in data_keys)

                        with tf.Session() as sess:
                            saver_finetune = tf.train.import_meta_graph(path_curr_net + "/" + file_name + ".meta")
                            saver_finetune.restore(sess, path_curr_net + "/" + file_name)
                        
                            # Trained variables of the first convolution layer
                            conv1_weights = {key : [v for v in tf.trainable_variables() if v.name == "conv1/weights_" + key + ":0"][0] for key in data_keys}
                            conv1_biases = {key : [v for v in tf.trainable_variables() if v.name == "conv1/biases_" + key + ":0"][0] for key in data_keys}
                            
                            # Extract results from ReLU layer
                            def relu_motif(data):
                                relu = {}
                                for i, key in enumerate(data_keys):
                                    conv1 = tf.nn.conv2d(data[i], 
                                                         conv1_weights[key],
                                                         strides = [1, 1, 1, 1],
                                                         padding = 'VALID')
                                    relu[key] = tf.nn.relu(tf.nn.bias_add(conv1, conv1_biases[key]))
                                return relu

                            relu_array = relu_motif(motif_data_node)
                            relu_array_background = relu_motif(motif_data_background_node)

                            # Motif signal of each motif detector on the dataset which are used for motif identification
                            motif_signal = sess.run([relu_array], feed_dict = {motif_data_node : motif_data})[0]
                            motif_signal_background = sess.run([relu_array_background], feed_dict = {motif_data_background_node: background_data})[0]
       
                    tf.reset_default_graph()    
                    
                    # Identify motifs for each feature type
                    for key in data_keys:
                        feature_dir = path_curr_net + "/" + key                    
                        util.make_dir(feature_dir)                        
                        os.chdir(feature_dir)
                        
                        f_binomial = open("Binomial.csv", 'w')
                        f_pvalue = open("Pvalue.csv", 'w')
                        
                        pvalue = np.zeros(shape = (FILTER_NUM, len(motif_cutoff)))
                        num_motif_instances = np.zeros(shape = (FILTER_NUM, len(motif_cutoff)))
                        mu_value = np.zeros(shape = (FILTER_NUM, len(motif_cutoff)))

                        # Calculate P-value
                        for i in range(motif_signal[key].shape[3]):
                            for j in range(len(motif_cutoff)):                               
                                mu = np.sum((np.max(motif_signal_background[key][..., i], axis = 1) > motif_cutoff[j] * np.max(motif_signal[key][..., i])) * 1) / (NUM_BACK_SEQ / motif_signal[key].shape[0])
                                x = np.sum((np.max(motif_signal[key][..., i], axis = 1) > motif_cutoff[j] * np.max(motif_signal[key][..., i])) * 1)

                                if mu == 0 and j != 0:
                                    pvalue[i, j] = np.inf
                                else:                                    
                                    p = mu / motif_signal[key].shape[0]
                                    if p > 1.0:
                                        p = 1.0
                                    pvalue[i, j] = binom_test(x, motif_signal[key].shape[0], p, alternative = "greater")
                                    
                                mu_value[i, j] = mu
                                num_motif_instances[i, j] = x                            
                                f_binomial.write(str(x) + '_' + str(mu) + '\t')
                                f_pvalue.write(str(pvalue[i, j]) + '\t')
                            f_binomial.write('\n')
                            f_pvalue.write('\n')
                        f_pvalue.close()
                        f_binomial.close()
                    
                        # Index of the minimum P-value for each motif detector across predefined motif cutoffs
                        index_min_pval = np.argmin(pvalue, axis = 1)

                        # Record the P-value of each motif
                        min_pvalue = np.amin(pvalue, axis = 1)
                        with open("min_motif_pvalue.csv", "w") as f:
                            f.write('Motif' + '\t' + 'Cutoff' + '\t' + "Max Motif Signal" + "\t" + 'P-value' + '\t' + '#Positive Instances' + "\t" + "#Negative Instances" + '\n')
                            for i in range(len(min_pvalue)):
                                f.write(str(i + 1) + '\t' + str(motif_cutoff[index_min_pval[i]]) + '\t' + str(np.max(motif_signal[key][..., i])) + "\t" + str(min_pvalue[i]) + '\t' + str(num_motif_instances[i, index_min_pval[i]]) + '\t' + str(mu_value[i, index_min_pval[i]]) + '\n')                                                
                        
                        # Align the sequence segments that maximally activate them 
                        for i in range(motif_signal[key].shape[3]):
                            motif_dir = feature_dir + "/motif_" + str(i + 1)                    
                            util.make_dir(motif_dir)                        
                            os.chdir(motif_dir)
                            
                            # Record the sequence segments in plain format and FASTA format
                            f_motif = open('seq_motif_instances.txt', 'w')
                            f_motifFA = open('seq_motif_instances.fa', 'w')
                            posi_motif = []    # Position of motif instances
                            
                            if feature_format == 'DNAShape':
                                shape_motif = []
                            
                            # The max activation in given sequences used for motif identification
                            max_activation = motif_cutoff[index_min_pval[i]] * np.max(motif_signal[key][..., i])
                            
                            # If motif signals are more enriched in positive sequences
                            if num_motif_instances[i, index_min_pval[i]] > mu_value[i, index_min_pval[i]]:
                                for j in range(motif_signal[key].shape[0]):
                                    if np.amax(motif_signal[key][j, ..., i]) > max_activation:
                                        index_max = np.argmax(motif_signal[key][j, ..., i])                                    
                                        f_motif.write(motif_sequences_alph[j][index_max : index_max + FILTER_LENGTH] + "\n")
                                        f_motifFA.write('>seq' + str(j + 1) + '\n')
                                        f_motifFA.write(motif_sequences_alph[j][index_max : index_max + FILTER_LENGTH] + "\n")
                                        posi_motif.append(index_max)

                                        if feature_format == 'DNAShape':
                                            shape_motif.append(util.convert_DNAShape(motif_data_comp[key][j, index_max : index_max + FILTER_LENGTH, 0, 0], feature_format, key))
                            
                            f_motifFA.close()
                            f_motif.close()
                            np.save("posi_motif.npy", np.array(posi_motif))
                            
                            if feature_format == 'DNAShape':
                                np.save("shape_motif_instances.npy", np.array(shape_motif))
                                np.savetxt("shape_motif_instances.csv", np.array(shape_motif), delimiter = ',')                            
                                
                            # Ignore empty seq_motif_instances
                            if os.stat('seq_motif_instances.fa').st_size != 0:                                                        
                                # Generate motif logo (normal and reverse complement) using weblogo 2.8.2
                                os.system(PATH_SEQLOGO + " -F PNG -a -n -Y -k 1 -c -w 50 -h 10 -f seq_motif_instances.fa > seqLogo.png")
                                os.system(PATH_SEQLOGO + " -F PNG -a -k 1 -c -w 50 -h 10 -f seq_motif_instances.fa > seqLogo_pure.png")

                        # Remove invalid motifs
                        os.chdir(feature_dir)
                        with open("min_motif_pvalue.csv") as f:
                            min_motif_pvalue = list(csv.reader(f, delimiter = "\t"))

                        for i in range(motif_signal[key].shape[3]):
                            flag = float(min_motif_pvalue[i + 1][4]) > float(min_motif_pvalue[i + 1][5]) and float(min_motif_pvalue[i + 1][4]) > 5 and float(min_motif_pvalue[i + 1][3]) < 0.0001
                            if not flag:
                                os.system("rm -rf " + feature_dir + "/motif_" + str(i + 1))
                            elif flag and feature_format == "Seq":
                                information_content = util.cal_IC("motif_" + str(i + 1) + "/seq_motif_instances.fa")
                                if not sum((np.array(information_content) > 1) * 1) >= 3:
                                    os.system("rm -rf " + feature_dir + "/motif_" + str(i + 1))
                            elif flag and feature_format == 'DNAShape':
                                shape_motif_instances = np.load("motif_" + str(i + 1) + "/shape_motif_instances.npy")
                                if not (np.any(np.mean(shape_motif_instances, axis=0) < SHAPE_THRESHOLD[KEY_SHAPE[key]][0]) or np.any(np.mean(shape_motif_instances, axis=0) > SHAPE_THRESHOLD[KEY_SHAPE[key]][1])):
                                    os.system("rm -rf " + feature_dir + "/motif_" + str(i + 1))
                                      
