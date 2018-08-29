'''
Train and test DESSO in predicting TF-DNA binding specifiticy
based on the 690 ENCODE ChIP-seq datasets

The 690 ChIP-seq datasets:
https://genome.ucsc.edu/ENCODE/downloads.html

TensorFlow install instructions:
https://tensorflow.org/get_started/os_setup.html
'''

from __future__ import division
import os
import sys

import numpy as np
import tensorflow as tf

from libs.constants import *
from libs.training import *
import libs.util as util

# Model parameters
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_integer('start_index', 0, 'Start index of the 690 ENCODE ChIP-seq datasets.')
flags.DEFINE_integer('end_index', 1, 'END index of the 690 ENCODE ChIP-seq datasets.')
flags.DEFINE_integer('peak_flank', 50, 'Number of flanking base pairs at each side of peak summit.')
flags.DEFINE_string('network', 'CNN', 'Neural network used in model training.')
flags.DEFINE_string('feature_format', 'Seq', 'Feature format of the input.')    # ['Seq', 'DNAShape']

PEAK_FLANK = [int(FLAGS.peak_flank)]
NETWORK = [FLAGS.network]
FEATURE_FORMAT = [FLAGS.feature_format]

# Flag for generating DNA shape features
if FLAGS.feature_format == 'Seq':
    DNA_SHAPE = False
elif FLAGS.feature_format == 'DNAShape':
    DNA_SHAPE = True
    
def main(argv = None):
    util.make_dir(PATH_OUTPUT)
    
    for peak_flank in PEAK_FLANK:
        print("Sequence length: %d" % (2 * peak_flank + 1))
        PATH_ENCODE = os.path.join(PATH_DATA, 'encode_' + str(2 * peak_flank + 1))       
        util.make_dir(os.path.join(PATH_OUTPUT, 'encode_' + str(2 * peak_flank + 1)))    
        
        for back_grou in BACK_GROU:
            PATH_OUTPUT_ENCODE = os.path.join(PATH_OUTPUT, 'encode_' + str(2 * peak_flank + 1), back_grou)
            util.make_dir(PATH_OUTPUT_ENCODE)
    
            # The name list of all ChIP-seq datasets
            [train_data_list, table_name_list] = util.get_data_name()

            for train_data_name in train_data_list[int(FLAGS.start_index) : int(FLAGS.end_index)]:
                print("Dataset: %s" % (train_data_name))
                table_name = table_name_list[train_data_list.index(train_data_name)]
                [peak_coor, peak_num] = util.get_peak_coor(PATH_ENCODE_TFBS_UNIF + "/" + table_name + ".narrowPeak.gz")

                # Record the results of the current training data
                path_curr_data = PATH_OUTPUT_ENCODE + '/' + train_data_name[-16:]
                util.make_dir(path_curr_data)     

                ############ Load Training Data ##############
                [train_sequences, train_targets, train_sequences_alph,  train_HelT, train_MGW, train_ProT, train_Roll, train_dnaShape] = \
                util.load_data_encode(PATH_ENCODE + "/" + train_data_name + "_AC.seq.gz", peak_coor, train_data_name, path_curr_data, "train", DNA_SHAPE, peak_flank, back_grou)
                
                ####### Convert data to one-hot format #######
                train_sequences_array = util.oneHot_data_encode(train_sequences)
                train_data_array = {}
                train_data_array['Seq'] = util.add_appendix(train_sequences_array, 0.25, FILTER_LENGTH)
                if DNA_SHAPE:
                    train_data_array['HelT'] = util.add_appendix(train_HelT[..., None], 0, FILTER_LENGTH)
                    train_data_array['MGW'] = util.add_appendix(train_MGW[..., None], 0, FILTER_LENGTH)  
                    train_data_array['ProT'] = util.add_appendix(train_ProT[..., None], 0, FILTER_LENGTH)
                    train_data_array['Roll'] = util.add_appendix(train_Roll[..., None], 0, FILTER_LENGTH)

                ############## Load Test Data ################
                [test_sequences, test_targets, test_sequences_alph, test_HelT, test_MGW, test_ProT, test_Roll, test_dnaShape] = \
                util.load_data_encode(PATH_ENCODE + "/" + train_data_name + "_B.seq.gz", peak_coor, train_data_name, path_curr_data, "test", DNA_SHAPE, peak_flank, back_grou)
                
                ####### Convert data to one-hot format #######
                test_sequences_array = util.oneHot_data_encode(test_sequences)
                test_data_array = {}
                test_data_array['Seq'] = util.add_appendix(test_sequences_array, 0.25, FILTER_LENGTH)
                if DNA_SHAPE:
                    test_data_array['HelT'] = util.add_appendix(test_HelT[..., None], 0, FILTER_LENGTH)
                    test_data_array['MGW'] = util.add_appendix(test_MGW[..., None], 0, FILTER_LENGTH)
                    test_data_array['ProT'] = util.add_appendix(test_ProT[..., None], 0, FILTER_LENGTH)
                    test_data_array['Roll'] = util.add_appendix(test_Roll[..., None], 0, FILTER_LENGTH)
                
                # Feature format
                for feature_format in FEATURE_FORMAT:
                    print('Feature format: %s' % (feature_format))
                    util.make_dir(path_curr_data + "/" + feature_format)
                    train_data_comp = util.prep_data_encode(train_data_array, feature_format)
                    test_data_comp = util.prep_data_encode(test_data_array, feature_format)
                    
                    # Neural network architecture
                    for net in NETWORK:
                        print('Network architecture: %s' % (net))
                        path_curr_net = path_curr_data + "/" + feature_format + "/" + net
                        util.make_dir(path_curr_net)

                        # Record hyperparameters and corresponding cross validation performance
                        samp_perf = np.zeros(shape = (NUM_SAMPLING, 6), dtype = np.float32)
                        print("***********************")
                        print("***** Calibration *****")
                        print("***********************")
                        
                        for sampling_count in range(NUM_SAMPLING):
                            print('Sampling count: %d' % (sampling_count))
                            # Select parameters from parameter space randomly                               
                            drop_out_rate = DROP_OUT_RATE[np.random.randint(len(DROP_OUT_RATE))]
                            learning_rate_init = np.random.uniform(5e-4, 0.05, 1)[0]
                            learning_momentum = np.random.uniform(0.95, 0.99, 1)[0]
                            lambda_regulate = np.random.uniform(1e-10, 1e-3, 1)[0]
                            std_weights = np.random.uniform(1e-5, 1e-1, 1)[0]

                            random_index = np.random.permutation(train_data_comp['0'].shape[0])
                            split_index = np.split(random_index, [train_data_comp['0'].shape[0] // FOLD_CV, 2 * (train_data_comp['0'].shape[0] // FOLD_CV), train_data_comp['0'].shape[0]])
                            cv_perf = np.zeros(shape = (3), dtype = np.float32)    # Performance of each cross validation

                            for cros_vali in range(FOLD_CV):
                                # Add appendix according to the fileter length
                                train_data = {key : np.delete(value, split_index[cros_vali], 0) for key, value in train_data_comp.items()} 
                                train_labels = np.delete(train_targets, split_index[cros_vali], 0)

                                validation_data = {key : value[split_index[cros_vali]] for key, value in train_data_comp.items()} 
                                validation_labels = train_targets[split_index[cros_vali]]

                                # Train a model based on given training data, validation data, and hyperparameter set
                                cv_perf[cros_vali] = train(train_data, train_labels, validation_data, validation_labels, path_curr_net, net, drop_out_rate, learning_rate_init, learning_momentum, lambda_regulate, std_weights)

                            samp_perf[sampling_count, :] = [drop_out_rate, learning_rate_init, learning_momentum, lambda_regulate, std_weights, np.mean(cv_perf)] 
                        
                        print("***********************")
                        print("******* Training ******")
                        print("***********************")
                        
                        # Permutate the training data and the corresponding targets
                        random_index = np.random.permutation(train_data_comp['0'].shape[0])
                        train_data_random = {key : value[random_index, ...] for key, value in train_data_comp.items()}
                        train_targets_random = train_targets[random_index, ...]

                        # 10% of the training data will be used as validation data
                        validation_size = train_data_comp['0'].shape[0] // 10
                        train_data = {key : value[validation_size:, ...] for key, value in train_data_random.items()}
                        train_labels = train_targets_random[validation_size:, ...]

                        validation_data = {key : value[:validation_size, ...] for key, value in train_data_random.items()}
                        validation_labels = train_targets_random[:validation_size, ...]

                        # Index of the best hyperparameter set
                        index_best_samp = np.argmax(samp_perf[:, 5])
                        _ = train(train_data, train_labels, validation_data, validation_labels, path_curr_net, net, \
                                   samp_perf[index_best_samp, 0], samp_perf[index_best_samp, 1], samp_perf[index_best_samp, 2], samp_perf[index_best_samp, 3], samp_perf[index_best_samp, 4], test_data_comp, test_targets)
                                
if __name__ == '__main__':
    tf.app.run()
