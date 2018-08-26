# DESSO (DEep Sequence and Shape mOtif)
DESSO is a deep learning-based framework that can be used to accurately identify both sequence and shape regulatory motifs from ENCODE ChIP-seq datasets.

<p align="center"> 
<img src="https://github.com/viyjy/DESSO/blob/master/Figure.PNG">
</p>

## Prerequisites and Dependencies
* Tensorflow 1.1.0
* CUDA 8.0.44
* Biopython 1.7.0
* Scikit-learn
* Download and 

## Model Training
Train CNN models on specified ChIP-seq datasets:
```
python DESSO.py --start_index 0 --end_index 1 --peak_flank 50 --network CNN --feature_format Seq
```
Arguments | Description
--------------|---------------------------------------------------------
--start_index | Start index of the 690 ENCODE ChIP-seq datasets
--end_index | END index of the 690 ENCODE ChIP-seq datasets
--peak_flank | Number of flanking base pairs at each side of peak summit
--network | Neural network used in model training
--feature_format | Feature format of the input

## Motif Prediction
Obtain either sequence or shape motif using the binomial distribution strategy based on the trained models above:

## Citation
If you use DESSO in your research, please cite the following paper.
