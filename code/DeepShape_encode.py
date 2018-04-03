from __future__ import division
import os
import sys
import numpy as np
from libs.encode_cons import *
from libs.global_cons import *
from libs.training import *
import libs.util as util
import tensorflow as tf

#########################
##### Main workflow #####
#########################
def main(argv = None):
    ### Hyperparameters initialization ###  
    drop_out_rate_list = [0.5, 0.75, 1]       # Dropout rate in hidden layer
        
    util.make_dir(PATH_OUTPUT)                # Make dir for output
    # PEAK_FLANK contains different integer values, each one indicates the size of flanking regions of peak summit
    for peak_flank in PEAK_FLANK:
        print("Sequence length: %d" % (2 * peak_flank + 1))
        PATH_ENCODE = os.path.join(PATH_DATA, 'encode_' + str(2 * peak_flank + 1))       # Path of peaks
        util.make_dir(os.path.join(PATH_OUTPUT, 'encode_' + str(2 * peak_flank + 1)))    # Path of output
        
        # BACK_GROU represents different background sequence generation methods
        for back_grou in BACK_GROU:
            print("Background sequence: %s" % (back_grou))
            PATH_OUTPUT_ENCODE = os.path.join(PATH_OUTPUT, 'encode_' + str(2 * peak_flank + 1), back_grou)
            util.make_dir(PATH_OUTPUT_ENCODE)
    
            # The name list of all ChIP-seq datasets
            [train_data_list, table_name_list] = util.get_data_name()

            for train_data_name in train_data_list[int(sys.argv[1]) : int(sys.argv[2])]:
                print("Dataset: %s" % (train_data_name))
                table_name = table_name_list[train_data_list.index(train_data_name)]
                peak_coor = util.get_peak_coor(PATH_ENCODE_TFBS_UNIF + "/" + table_name + ".narrowPeak.gz", peak_flank)

                # Record the results of the current training data
                path_curr_data = PATH_OUTPUT_ENCODE + '/' + train_data_name[-16:]
                util.make_dir(path_curr_data)     

                ################# Load Data ##################
                # All datasets are loaded using the following function
                [train_sequences, train_targets, train_sequences_alph,  train_HelT, train_MGW, train_ProT, train_Roll, train_dnaShape] = \
                util.load_data_encode(PATH_ENCODE + "/" + train_data_name + "_AC.seq.gz", peak_coor, train_data_name, path_curr_data, "train", DNA_SHAPE, peak_flank, back_grou)
                
                ####### Convert data to one-hot format #######
                [train_sequences_array, train_targets_array, train_HelT_array, train_MGW_array, train_ProT_array, train_Roll_array] = \
                util.oneHot_data_encode(train_sequences, train_targets, train_HelT, train_MGW, train_ProT, train_Roll, DNA_SHAPE)
                train_data_array = {}
                train_data_array['Seq'] = util.add_appendix(train_sequences_array, 0.25, filter_length)
                train_data_array['HelT'] = util.add_appendix(train_HelT_array, 0, filter_length)
                train_data_array['MGW'] = util.add_appendix(train_MGW_array, 0, filter_length)  
                train_data_array['ProT'] = util.add_appendix(train_ProT_array, 0, filter_length)
                train_data_array['Roll'] = util.add_appendix(train_Roll_array, 0, filter_length)
                
                # All the data are loaded using the following function.
                [test_sequences, test_targets, test_sequences_alph, test_HelT, test_MGW, test_ProT, test_Roll, test_dnaShape] = \
                util.load_data_encode(PATH_ENCODE + "/" + train_data_name + "_B.seq.gz", peak_coor, train_data_name, path_curr_data, "test", DNA_SHAPE, peak_flank, back_grou)
                
                ####### Convert data to one-hot format #######
                [test_sequences_array, test_targets_array, test_HelT_array, test_MGW_array, test_ProT_array, test_Roll_array] = \
                util.oneHot_data_encode(test_sequences, test_targets, test_HelT, test_MGW, test_ProT, test_Roll, DNA_SHAPE)
                test_targets = test_targets_array
                test_data_array = {}
                test_data_array['Seq'] = util.add_appendix(test_sequences_array, 0.25, filter_length)
                test_data_array['HelT'] = util.add_appendix(test_HelT_array, 0, filter_length)
                test_data_array['MGW'] = util.add_appendix(test_MGW_array, 0, filter_length)
                test_data_array['ProT'] = util.add_appendix(test_ProT_array, 0, filter_length)
                test_data_array['Roll'] = util.add_appendix(test_Roll_array, 0, filter_length)
                
                # Different feature format
                for feature_format in FEATURE_FORMAT:
                    print('Feature format: %s' % (feature_format))
                    util.make_dir(path_curr_data + "/" + feature_format)
                    train_data_comp = util.prep_data_encode(train_data_array, feature_format)
                    test_data_comp = util.prep_data_encode(test_data_array, feature_format)
                    
                    # Different neural network architecture
                    for net in NETWORK:
                        print('Network architecture: %s' % (net))
                        path_curr_net = path_curr_data + "/" + feature_format + "/" + net
                        util.make_dir(path_curr_net)

                        # Record hyperparameters and corresponding cross validation performance (mean)
                        samp_perf = np.zeros(shape = (NUM_SAMPLING, 6), dtype = np.float32)

                        # Iterate NUM_SAMPLING times using different hyperparameter sets in each iteration. 
                        # Each hyperparameter set will be used to train a model using three-fold cross validation strategy based on whole training data
                        # Select the best hyperparameter set based on the average performance of each three-fold cross validation
                        # The selected best hyperparameter set will be used to train a final model based on whole training data, the apply the trained model on test data
                        print("***********************")
                        print("***** Calibration *****")
                        print("***********************")
                        
                        for sampling_count in range(NUM_SAMPLING):
                            print('Sampling count: %d' % (sampling_count))
                            # Select parameters from parameter space randomly                               
                            drop_out_rate = drop_out_rate_list[np.random.randint(len(drop_out_rate_list))]
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
                                train_targets = np.delete(train_targets_array, split_index[cros_vali], 0)

                                validation_data = {key : value[split_index[cros_vali]] for key, value in train_data_comp.items()} 
                                validation_targets = train_targets_array[split_index[cros_vali]]

                                # Train a model based on given training data, validation data, and hyperparameter set
                                cv_perf[cros_vali] = train(train_data, train_targets, validation_data, validation_targets, path_curr_net, net, drop_out_rate, learning_rate_init, learning_momentum, lambda_regulate, std_weights)

                            samp_perf[sampling_count, :] = [drop_out_rate, learning_rate_init, learning_momentum, lambda_regulate, std_weights, np.mean(cv_perf)] 
                        
                        print("***********************")
                        print("******* Training ******")
                        print("***********************")
                        
                        # Permutate the training data and training targets
                        random_index = np.random.permutation(train_data_comp['0'].shape[0])
                        train_data_random = {key : value[random_index, ...] for key, value in train_data_comp.items()}
                        train_targets_random = train_targets_array[random_index, ...]

                        # 10% of the training data will be used as validation data
                        validation_size = train_data_comp['0'].shape[0] // 10
                        train_data = {key : value[validation_size:, ...] for key, value in train_data_random.items()}
                        train_targets = train_targets_random[validation_size:, ...]

                        validation_data = {key : value[:validation_size, ...] for key, value in train_data_random.items()}
                        validation_targets = train_targets_random[:validation_size, ...]

                        # Index of the best hyperparameter set
                        index_best_samp = np.argmax(samp_perf[:, 5])
                        _ = train(train_data, train_targets, validation_data, validation_targets, path_curr_net, net, \
                                   samp_perf[index_best_samp, 0], samp_perf[index_best_samp, 1], samp_perf[index_best_samp, 2], samp_perf[index_best_samp, 3], samp_perf[index_best_samp, 4], test_data_comp, test_targets)
                
if __name__ == '__main__':
    tf.app.run()
