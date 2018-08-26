from __future__ import division
import sys
import time

from constants import *
from model import *
import util as util

import numpy as np
import tensorflow as tf
sys.path.append("/home/xsede/users/xs-viyjy/lib/python2.7/site-packages/")
from sklearn.metrics import roc_auc_score

def train(train_data, train_targets, validation_data, validation_targets, path_curr_net, net, drop_out_rate, learning_rate_init, learning_momentum, lambda_regulate, std_weights, test_data = None, test_targets = None):
    best_validation_perf = 0    # Best performance of the trained model on validation data
    
    # Set the current graph as default at each time
    with tf.Graph().as_default():
        train_size = train_data['0'].shape[0]    # Size of the training data

        ### Placeholder inputs which will be fed into graph ###
        data_keys = train_data.keys()
        train_data_node = tuple(tf.placeholder(tf.float32, shape = (BATCH_SIZE, train_data[key].shape[1], train_data[key].shape[2], train_data[key].shape[3])) for key in data_keys)
        train_targets_node = tf.placeholder(tf.float32, shape = (BATCH_SIZE,))
        eval_data_node = tuple(tf.placeholder(tf.float32, shape = (EVAL_BATCH_SIZE, validation_data[key].shape[1], validation_data[key].shape[2], validation_data[key].shape[3])) for key in data_keys)

        ### Trainable parameters ###
        # These parameters will be optimized in the training process to minimize the loss function
        conv1_weights = {}
        conv1_biases = {}
        with tf.variable_scope('conv1', reuse = True):
            for key in data_keys:
                conv1_weights[key] = tf.Variable(tf.truncated_normal([FILTER_LENGTH, train_data[key].shape[2], train_data[key].shape[3], FILTER_NUM], stddev = std_weights, dtype = tf.float32, seed = SEED), name = ('weights_' + key))
                conv1_biases[key] = tf.Variable(tf.constant(0, shape = [FILTER_NUM], dtype = tf.float32), name = ('biases_' + key))
        
        fc1_weights, fc1_biases = util.para_trainable([len(data_keys) * FILTER_NUM, HIDDEN_LAYER_SIZE], std_weights, 0, 'fc1')
        fc2_weights, fc2_biases = util.para_trainable([HIDDEN_LAYER_SIZE, NUM_LABELS], std_weights, 0, 'fc2')
        para = {'conv1_weights': conv1_weights,
                'conv1_biases': conv1_biases,
                'fc1_weights': fc1_weights,
                'fc1_biases': fc1_biases,
                'fc2_weights': fc2_weights,
                'fc2_biases': fc2_biases,
                'drop_out_rate': drop_out_rate,
                'std_weights': std_weights}

        ### Optimization ###
        logits = model(train_data_node, para, data_keys, net, True)
        loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels = train_targets_node, logits = tf.reshape(logits, [-1])))                              

        # L2 regularization #
        regularizers = (tf.nn.l2_loss(fc1_weights) + tf.nn.l2_loss(fc1_biases) + tf.nn.l2_loss(fc2_weights) + tf.nn.l2_loss(fc2_biases))

        # Loss function
        loss += (lambda_regulate * regularizers)

        ### Optimizer ###
        batch = tf.Variable(0, dtype = tf.float32)
        learning_rate = tf.train.exponential_decay(learning_rate_init, batch * BATCH_SIZE, train_size, 0.95, staircase = True)

        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
        with tf.control_dependencies(update_ops):
            optimizer = tf.train.MomentumOptimizer(learning_rate, learning_momentum).minimize(loss, global_step = batch)

        # Prediction of each minibatch
        train_predictions = tf.nn.sigmoid(logits)

        # Prediction of validation data
        eval_prediction = tf.nn.sigmoid(model(eval_data_node, para, data_keys, net))

        def eval_in_batches(data, sess):
            size = data['0'].shape[0]
            if size < EVAL_BATCH_SIZE:
                raise ValueError("Batch size larger than dataset: %d" % size)
            predictions = np.ndarray(shape = (size, NUM_LABELS), dtype = np.float32)
            for begin in xrange(0, size, EVAL_BATCH_SIZE):
                end = begin + EVAL_BATCH_SIZE
                if end <= size:
                    predictions[begin:end, :] = sess.run(
                        eval_prediction,
                        feed_dict = {eval_data_node : tuple(data[key][begin : end, ...] for key in data_keys)})
                else:
                    batch_predictions = sess.run(
                        eval_prediction,
                        feed_dict = {eval_data_node: tuple(data[key][-EVAL_BATCH_SIZE:, ...] for key in data_keys)})
                    predictions[begin:, :] = batch_predictions[begin - size:, :]
            return predictions

        #####################################################
        ########## Local session to train the model #########
        saver = tf.train.Saver()    # Save the model with maximum AUC
        start_time = time.time()    # Record the time for each update
        with tf.Session() as sess:
            tf.global_variables_initializer().run()

            # Early-stopping strategy based on a geometrically increasing amount of patience
            eval_frequency = train_size // BATCH_SIZE     # Evaluate the trained model after each epoch
            patience = 5000                               
            patience_increase = 2 
            iteration_count = 0                           # Count the parameter update time
            epoch = 0
            done_looping = False

            while epoch < NUM_EPOCHS and (not done_looping):
                epoch += 1
                random_index = np.random.permutation(train_size)    # Permute the batch_data and batch_targets cross each epoch
                for start_point in xrange(0, train_size - BATCH_SIZE, BATCH_SIZE):
                    # Extract the batch data and targets into feed_dict which will be fed into the graph 
                    batch_data = tuple(train_data[key][random_index[start_point : start_point + BATCH_SIZE], ...] for key in data_keys)
                    batch_targets = train_targets[random_index[start_point : start_point + BATCH_SIZE], ...]
                    feed_dict = {train_data_node : batch_data, train_targets_node : batch_targets}

                    _, curr_loss, curr_learning_rate, predictions = sess.run([optimizer, loss, learning_rate, train_predictions], feed_dict = feed_dict)

                    ### Evaluate the current model on validation data 
                    if (iteration_count + 1) % eval_frequency == 0:
                        elapsed_time = time.time() - start_time
                        start_time = time.time()											                                                        

                        ### Evaluate the performance of current model on validation data underlying the AUC 						
                        validation_AUC = roc_auc_score(validation_targets, eval_in_batches(validation_data, sess).flatten())
                        print('Epoch: %d, Iteration: %d, Time: %d seconds, Minibatch validation AUC: %f' % (epoch, iteration_count + 1, elapsed_time, validation_AUC))						

                        ### If current AUC is larger than best_validation_perf
                        if validation_AUC > best_validation_perf:
                            best_validation_perf = validation_AUC    # update the max AUC
                            patience = max(patience, iteration_count * patience_increase)

                            if test_data is not None:
                                test_AUC = roc_auc_score(test_targets, eval_in_batches(test_data, sess).flatten())
                                with open(path_curr_net + "/Test_result.txt", 'w') as f:
                                    f.write(str(test_AUC))
                                # Save the trained model                                                                        
                                save_path = saver.save(sess, path_curr_net + "/" + "model.ckpt")
                        sys.stdout.flush()

                    if patience <= iteration_count:													
                        done_looping = True
                        break
                    iteration_count += 1
    tf.reset_default_graph()
    return best_validation_perf
