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
* Download [GRCh37.p13.genome.fa](http://bmbl.sdstate.edu/downloadFiles/GRCh37.p13.genome.fa) and [encode_101_background](http://bmbl.sdstate.edu/downloadFiles/encode_101_background.7z), then put them into ```data/```

## Model Training
Train CNN models for specificied datasets: 
```
cd code/
python DESSO.py --start_index 0 --end_index 1 --peak_flank 50 --network CNN --feature_format Seq
```
Arguments | Description
--------------|---------------------------------------------------------
--start_index | Start index of the 690 ENCODE ChIP-seq datasets
--end_index | END index of the 690 ENCODE ChIP-seq datasets
--peak_flank | Number of flanking base pairs at each side of peak summit
--network | Neural network used in model training
--feature_format | Feature format of the input

DESSO can be applied to the [690 ChIP-seq datasets](https://genome.ucsc.edu/ENCODE/downloads.html). 
```--start_index 0 --end_index 1``` above indicates the first dataset (i.e., wgEncodeEH002288) was used for model training.
```--peak_flank 50``` indicates the peak length is 101 bps

## Motif Prediction
Obtain either sequence or shape motif using the binomial distribution strategy based on the trained models above:

## Citation
If you use DESSO in your research, please cite the following paper.
