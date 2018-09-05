# DESSO (DEep Sequence and Shape mOtif)
DESSO is a deep learning-based framework that can be used to accurately identify both sequence and shape regulatory motifs from the human genome. The performance of DESSO was evaluated on the [690 ChIP-seq datasets](https://genome.ucsc.edu/ENCODE/downloads.html).

<p align="center"> 
<img src="https://github.com/viyjy/DESSO/blob/master/workflow.PNG">
</p>

## Prerequisites and Dependencies
* Tensorflow 1.1.0 [[Install]](https://www.tensorflow.org/install/)
* CUDA 8.0.44
* Python 2.7
* Biopython 1.7.0
* Scikit-learn
* Download [GRCh37.p13.genome.fa](http://bmbl.sdstate.edu/DESSO/tools/GRCh37.p13.genome.fa.zip) and [encode_101_background](http://bmbl.sdstate.edu/DESSO/tools/encode_101_background.zip), then unzip them and put them into ```data/```
* ```data/encode_101```, ```data/encode_1001```, and ```data/TfbsUniform_hg19_ENCODE``` only contain wgEncodeEH002288-related data as an example. The source code and whole data can be accessed at [code+whole data](http://bmbl.sdstate.edu/DESSO/tools/DESSO-master-whole.zip).

## Model Training Based on CNN
Train CNN models on specified datasets: 
```
cd code/
python train.py --start_index 0 --end_index 1 --peak_flank 50 --network CNN --feature_format Seq
```
Arguments | Description
--------------|---------------------------------------------------------
--start_index | Start index of the 690 ENCODE ChIP-seq datasets
--end_index | END index of the 690 ENCODE ChIP-seq datasets
--peak_flank | Number of flanking base pairs at each side of peak summit (default is 50)
--network | Neural network used in model training (default is CNN)
--feature_format | Feature format of the input (default is Seq)

```--start_index 0 --end_index 1``` indicates the first dataset (i.e., wgEncodeEH002288). E.g., to train models for the second and third datasets, use ```--start_index 1 --end_index 3``` <br/>
```--peak_flank 50``` indicates the peak length is (2 * 50 + 1) = 101 base pairs <br/>
```--network``` indicates that CNN is used here <br/>
```--feature_format``` can be Seq or DNAShape, where Seq indicates the input is DNA sequences, DNAShape indicates the input is the combination of four DNA shape features (i.e., HelT, MGW, ProT, and Roll).

## Motif Prediction
Obtain either sequence or shape motifs based on the trained models above:
```
cd code/
python predict.py --start_index 0 --end_index 1 --peak_flank 50 --network CNN --feature_format Seq --start_cutoff 0.01 --end_cutoff 1 --step_cutoff 0.03
```
Arguments | Description
----------|----------------------------------------------------------
--start_cutoff | Start of the motif cutoff interval (default is 0.01)
--end_cutoff | End of the motif cutoff interval (default is 1)
--step_cutoff | Increament of the cutoff (default is 0.03)

```--feature_format Seq``` indicates that sequence motifs will be predicted. To identify shape motifs, use ```--feature_format DNAShape``` instead.

## Predict TF-DNA Binding Specicitity Using Gated-CNN and Long DNA Sequence
```
cd code/
python train.py --start_index 0 --end_index 1 --peak_flank 500 --network GCNN --feature_format Seq
```
```--peak_flank 500``` indicates the peak length is (2 * 500 + 1) = 1001 base pairs <br/>

## Citation
If you use DESSO in your research, please cite the following paper.
