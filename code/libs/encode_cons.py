### The global constants used in related script ###
from global_cons import *
'''
PEAK_FLANK = [50, 500]                    # Number of flanking bps at each side of peak summit
BACK_GROU = ['rand_geno', 'dinu_shuf']    # Methods for generating background (or negative) sequences:
                                          # 'dinu_shuf' indicates dinucleotide-preserving shuffle,
                                          # 'rand_geno' indicates random genome segments selection
NETWORK = ['CNN', 'GCNN']                 # Type of network: CNN (convolutional neural network) and GCNN (Gated CNN)
'''
PEAK_FLANK = [50]
BACK_GROU = ['rand_geno']
NETWORK = ['CNN']
FEATURE_FORMAT = ['Seq', 'HelT', 'MGW', 'ProT', 'Roll', 'DNAShape', 'Seq_DNAShape']
#FEATURE_FORMAT = ['Seq']

DNASHAPE_FEATURE = ['HelT', 'MGW', 'ProT', 'Roll']
#PATH_NBT_S7 = PATH_DATA + '/deepbind/nbt.3300-S7.txt'
PATH_ENCODE_TFBS = PATH_DATA + "/encode_tfbs.txt"
PATH_ENCODE_TFBS_UNIF = PATH_DATA + "/TfbsUniform_hg19_ENCODE"

# Different combination of the available features, including 'S' indicates sequences,
# 'M' indicates MGW, 'P' indicates ProT, 'SM' indicates sequneces and MGW,
# 'SP' indicates sequences and ProT, 'MP' indicates MGW and ProT, 'SMP' indicates sequences, MGW and ProT
#FEATURE_FORMAT = ['S', 'M', 'P', 'SM', 'SP', 'MP', 'SMP']

num_channels = 1
filter_length = 24        # The lengh of each filter in convolutional layer 
filter_num = 16           # The number of the filters in convolutional layer 
hidden_layer_size = 32    # The number of neurons in hidden layer
NUM_SHUFFLE = 1000         # shuffle time used in motif cutoff determination
SEQUENCE_WIDTH = 4        # 'ACGT'

#ORIGINAL_SEQUENCE_LENGTH = 2 * PEAK_FLANK + 1    # The original sequence length
#VALID_SEQUENCE_LENGTH = 2 * PEAK_FLANK + 1       # The valid sequence length since the last five bases should be eliminated
#VALIDATION_SIZE = 2000                           # Size of the validation data set
FOLD_CV = 3                                       # The fold of cross validation
SEED = 66478                # The seed of random number generator
NUM_LABELS = 1              # Binary classification
NUM_SAMPLING = 10            # Hyperparameters sampling time
NUM_EPOCHS = 30             # Epoches for training
BATCH_SIZE = 64             # Batch size in data training
EVAL_BATCH_SIZE = 64        # Batch size in model evaluation
MINIMUM_POSI_SEQ = 10000    # Minimum positive sequences
DNA_SHAPE = True           # DNAShape usage flag

# Three human motif databases including JASPAR (386 Motifs), HOCOMOCO (641 Motifs), and TRANSFAC (208 Motifs)
JASPAR_DB = PATH_MOTIF_DB + "/JASPAR/JASPAR_CORE_2016_vertebrates.meme"
HOCOMOCO_DB = PATH_MOTIF_DB + "/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
TRANSFAC_DB = PATH_MOTIF_DB + "/TRANSFAC/TRANSFAC_2005_HUMAN.meme"
