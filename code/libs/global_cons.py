### The global constants ###
import os
PATH_PROJ1 = os.getcwd()[:-5]                            # The path of this project

# The main filders used in storing code, data and output
PATH_CODE = os.path.join(PATH_PROJ1, 'code')             # The path of code used in this project
PATH_DATA = os.path.join(PATH_PROJ1, 'data')             # The path of data used in this project
PATH_OUTPUT = os.path.join(PATH_PROJ1, 'output')         # The path of the obtained output 

# Folders included in code/
PATH_SVR = os.path.join(PATH_CODE, 'SVR')                # The path of SVR-related files
PATH_LIBS = os.path.join(PATH_CODE, 'libs')              # The path of libs used for storing tools and function

# Folders and executable files included in code/libs/
PATH_MOTIF_DB = os.path.join(PATH_LIBS, 'motif_databases')            # The path of the motif databases downloaded from MEME
PATH_SEQLOGO = os.path.join(PATH_LIBS, 'weblogo.2.8.2/seqlogo')           # The path of weblogo2 used for motif logo generation

#PATH_DEEPBIND = os.path.join(PATH_DATA, 'deepbind-v0.11-linux/deepbind')    # The path of the trained models of DeepBind
'''
PATH_OUTPUT_SVR = os.path.join(PATH_OUTPUT, 'SVR')    # The output path of SVR
PATH_OUTPUT_DEEPSHAPE = os.path.join(PATH_OUTPUT, 'deepShape')        # The output path of DeepShape
PATH_OUTPUT_DREAM5 = os.path.join(PATH_OUTPUT_DEEPSHAPE, 'dream5')    # The output path of DeepShape based on DREAM5
PATH_OUTPUT_ENCODE = os.path.join(PATH_OUTPUT_DEEPSHAPE, 'encode')    # The output path of DeepShape based on ENCODE
'''

# The base-num pairs
BASE_NUM = {'A' : 1, 
            'C' : 2,
            'G' : 3,
            'T' : 4,
            'N' : 0
           }
