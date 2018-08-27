### The global constants ###
import os

PATH_DESSO = os.getcwd()[:-5]                            
PATH_CODE = os.path.join(PATH_DESSO, 'code')             
PATH_DATA = os.path.join(PATH_DESSO, 'data')             
PATH_OUTPUT = os.path.join(PATH_DESSO, 'output')         
PATH_LIBS = os.path.join(PATH_CODE, 'libs')     

PATH_ENCODE_TFBS = PATH_DATA + "/encode_tfbs.txt"
PATH_ENCODE_TFBS_UNIF = PATH_DATA + "/TfbsUniform_hg19_ENCODE"
PATH_SEQLOGO = os.path.join(PATH_LIBS, 'weblogo.2.8.2/seqlogo')

BACK_GROU = ['gc_match']    # Methods for generating negative sequences ['gc_match', 'rand_geno', 'dinu_shuf']
SEQUENCE_WIDTH = 4          # 'ACGT'
NUM_BACK_SEQ = 500000       # Number of background sequences
MINIMUM_POSI_SEQ = 10000    # Minimum positive sequences

# Basic model parameters
FILTER_LENGTH = 24          # Lengh of each filter in convolutional layer 
FILTER_NUM = 16             # Number of filters in convolutional layer 
HIDDEN_LAYER_SIZE = 32      # Number of neurons in hidden layer
FOLD_CV = 3                 # Fold of cross validation
SEED = 66478                # Seed of random number generator
NUM_LABELS = 1              # Binary classification
NUM_SAMPLING = 10           # Hyperparameters sampling time
NUM_EPOCHS = 30             # Epoches for training
BATCH_SIZE = 64             # Batch size in model training
EVAL_BATCH_SIZE = 64        # Batch size in model evaluation
DROP_OUT_RATE = [0.5, 0.75, 1]    # Dropout rate in hidden layer

# The base-num pairs
BASE_NUM = {'A' : 1, 
            'C' : 2,
            'G' : 3,
            'T' : 4,
            'N' : 0
           }

DNASHAPE_FEATURE = ['HelT', 'MGW', 'ProT', 'Roll']

# These constants are used to scale DNA shape features
MGW_max_first = 6.2
MGW_min_first = 2.85
MGW_span_first = MGW_max_first - MGW_min_first

ProT_max_first = -0.03
ProT_min_first = -16.51
ProT_span_first = ProT_max_first - ProT_min_first

Roll_max_first = 8.64
Roll_min_first = -8.57
Roll_span_first = Roll_max_first - Roll_min_first

HelT_max_first = 38.05
HelT_min_first = 30.94
HelT_span_first = HelT_max_first - HelT_min_first
