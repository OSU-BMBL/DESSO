### The global constants used in other script
# The path of sequences, targets, probe_biases and TF_ID of PBM data in DREAM5 challenge
PATH_SEQUENCE = '../data/dream5/pbm/sequences.tsv.gz'
PATH_TARGETS = '../data/dream5/pbm/targets.tsv.gz'
PATH_PROBE_BIASES = '../data/dream5/pbm/probe_biases.npz'
PATH_TF_ID = '../data/dream5/pbm/tfids.txt'

# The path of DNA shape information of sequences in PBM
PATH_HelT = '../data/dream5/pbm/DNAShape/all/binding_sites.HelT'    # 
PATH_MGW = '../data/dream5/pbm/DNAShape/all/binding_sites.MGW'    # All the items in MGW are larger than and equal to zero
PATH_ProT = '../data/dream5/pbm/DNAShape/all/binding_sites.ProT'    # All the items in ProT are less than and equal to zero
PATH_Roll = '../data/dream5/pbm/DNAShape/all/binding_sites.Roll'    # 

# Different combination of the available features, including 'S' indicates sequences, 'M' indicates MGW, 'P' indicates ProT, 'SM' indicates sequneces and MGW, 'SP' indicates sequences and ProT, 'MP' indicates MGW and ProT, 'SMP' indicates sequences, MGW and ProT
#FEATURE_FORMAT = ['S', 'M', 'P', 'SM', 'SP', 'MP', 'SMP']
FEATURE_FORMAT = ['S']
SEQUENCE_WIDTH = 4    # 'ACGT'
ORIGINAL_SEQUENCE_LENGTH = 40    # The original sequence length
VALID_SEQUENCE_LENGTH = 35    # The valid sequence length since the last five bases should be eliminated

VALIDATION_SIZE = 5000     # Size of the validation data set
#SEED = 66478
NUM_SAMPLING = 30    # Hyperparameters sampling time
NUM_EPOCHS = 30    # Epoches for training
BATCH_SIZE = 64    # Batch size in data training
EVAL_BATCH_SIZE = 64    # Batch size in model evaluation
EVAL_FREQUENCY = 100    # Frequency of model evaluation
LAMBDA = 5e-4    # Regularized rate

