''' Models with different neural network architectures '''
from __future__ import division
from libs.constants import *
import libs.util as util
import tensorflow as tf

def model(data, para, data_keys, net, train = False):
    '''
    Args:
        data: Model input
        para: Trainable parameters
        data_keys: Key of the data, each key represents a feature format
        net: Neural network architecture, i.e., CNN or GCNN
        train: A flag representing model training or model test
    Return:
        Logits derived from the selected model
    '''
    
    ''' Gated Convolutional Neural Network '''
    if net == 'GCNN':        
        conv1 = tf.nn.conv2d(data, 
                            para['conv1_weights'],
                            strides = [1, 1, 1, 1],
                            padding = 'VALID')
        relu1 = tf.nn.relu(tf.nn.bias_add(conv1, para['conv1_biases']))
        pool1 = tf.nn.max_pool(relu1,
                               ksize = [1, 10, 1, 1],
                               strides = [1, 10, 1, 1],
                               padding = 'VALID')

        last_pool_shape = pool1.get_shape().as_list()
        motif_embed = tf.reshape(pool1, [last_pool_shape[0], last_pool_shape[1], last_pool_shape[3], last_pool_shape[2]])                
    
        num_layers = 3
        filter_size = FILTER_NUM
        filter_h = 5
        filter_w = last_pool_shape[3]
        block_size = 2

        h, res_input = motif_embed, motif_embed

        for i in range(num_layers):
            h = tf.pad(h, [[0, 0], [filter_h - 1, 0], [0, 0], [0, 0]], 'CONSTANT')
            fanin_depth = h.get_shape().as_list()[-1]
            shape = (filter_h, filter_w, fanin_depth, filter_size)

            with tf.variable_scope("layer_%d"%i) as scope:
                try:    
                    linear_W = tf.get_variable("linear_W", shape, tf.float32, tf.random_normal_initializer(0.0, para['std_weights'], seed=SEED))
                    linear_b = tf.get_variable("linear_b", shape[-1], tf.float32, tf.constant_initializer(0))

                    gated_W = tf.get_variable("gated_W", shape, tf.float32, tf.random_normal_initializer(0.0, para['std_weights'], seed=SEED))
                    gated_b = tf.get_variable("gated_b", shape[-1], tf.float32, tf.constant_initializer(0))

                except ValueError:
                    scope.reuse_variables()
                    linear_W = tf.get_variable("linear_W")
                    linear_b = tf.get_variable("linear_b")

                    gated_W = tf.get_variable("gated_W")
                    gated_b = tf.get_variable("gated_b")

                conv_w = tf.nn.bias_add(tf.nn.conv2d(h, linear_W, strides = [1, 1, 1, 1], padding = 'VALID'), linear_b)
                conv_v = tf.nn.bias_add(tf.nn.conv2d(h, gated_W, strides = [1, 1, 1, 1], padding = 'VALID'), gated_b)
                h = conv_w * tf.nn.sigmoid(conv_v)

                if i == (num_layers - 1):
                    h = tf.nn.max_pool(h,
                                       ksize = [1, h.get_shape().as_list()[1], 1, 1],
                                       strides = [1, h.get_shape().as_list()[1], 1, 1],
                                       padding = 'VALID')

                h_shape = h.get_shape().as_list()
                h = tf.reshape(h, [h_shape[0], h_shape[1], h_shape[3], h_shape[2]])
                if (i + 1) % block_size == 0:
                    h += res_input
                    res_input = h


        h_shape = h.get_shape().as_list()
        fc_input = tf.reshape(h, [h_shape[0], h_shape[1] * h_shape[2] * h_shape[3]])
        hidden1 = tf.nn.relu(tf.matmul(fc_input, para['fc1_weights']) + para['fc1_biases'])
        if train:
            hidden1 = tf.nn.dropout(hidden1, para['drop_out_rate'], seed=SEED)

        return tf.matmul(hidden1, para['fc2_weights']) + para['fc2_biases']
    elif net == 'CNN':
        ''' Convolutional neural network '''
        for i, key in enumerate(data_keys):
            conv1 = tf.nn.conv2d(data[i], 
                                 para['conv1_weights'][key],
                                 strides = [1, 1, 1, 1],
                                 padding = 'VALID')
            relu1 = tf.nn.relu(tf.nn.bias_add(conv1, para['conv1_biases'][key]))
            pool1 = tf.nn.max_pool(relu1,
                                   ksize = [1, data[i].get_shape().as_list()[1] - FILTER_LENGTH + 1, 1, 1],
                                   strides = [1, 1, 1, 1],
                                   padding = 'VALID')

            # Fully connected layer #
            last_pool_shape = pool1.get_shape().as_list()
            if i == 0:
                fc_input = tf.reshape(pool1, [last_pool_shape[0], last_pool_shape[1] * last_pool_shape[2] * last_pool_shape[3]])
            else:
                fc_input = tf.concat([fc_input, tf.reshape(pool1, [last_pool_shape[0], last_pool_shape[1] * last_pool_shape[2] * last_pool_shape[3]])], 1)

        hidden1 = tf.nn.relu(tf.matmul(fc_input, para['fc1_weights']) + para['fc1_biases'])
        if train:
            hidden1 = tf.nn.dropout(hidden1, para['drop_out_rate'], seed=SEED)

        return tf.matmul(hidden1, para['fc2_weights']) + para['fc2_biases']
