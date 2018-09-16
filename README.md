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
* ```data/encode_101```, ```data/encode_1001```, and ```data/TfbsUniform_hg19_ENCODE``` only contain wgEncodeEH002288-related data as an example, owing to the file size limit. To access the source code and whole datasets (totally about 5.9GB) without additional manipulation, just click on [code+whole data](http://bmbl.sdstate.edu/DESSO/tools/DESSO-master-whole.zip).

## Model Training Based on Convolutional Neural Network (CNN)
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

```--start_index 0 --end_index 1``` indicates the first dataset (i.e., wgEncodeEH002288). For example, to train models for the second and third datasets, use ```--start_index 1 --end_index 3``` <br/>
```--peak_flank 50``` indicates the peak length is (2 * 50 + 1) = 101 base pairs <br/>
```--network``` indicates that CNN is used here <br/>
```--feature_format``` can be Seq or DNAShape, where Seq indicates the input is DNA sequences, DNAShape indicates the input is the combination of four DNA shape features (i.e., HelT, MGW, ProT, and Roll).

### Output
If ```--feature_format Seq``` was used, the trained model can be found at ```/output/encode_101/gc_match/wgEncodeEH002288/Seq/CNN```, together with ```Test_result.txt``` indicating the area under the receiver operating characteristic curve (AUC) of the trained model in predicting TF-DNA binding specificity on the test data. <br/>
If ```--feature_format DNAShape``` was used, the trained model is located at ```/output/encode_101/gc_match/wgEncodeEH002288/DNAShape/CNN```.

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

### Output
For ```--feature_format Seq```, the predicted sequence motifs are in ```output/encode_101/gc_match/wgEncodeEH002288/Seq/CNN/0```. <br/>
For ```--feature_format DNAShape```, four kinds of shape motifs would be predicted as shown in the following table:

Location | Type of predicted shape motif
-----------------------------------------------------------------|-----------------------------
output/encode_101/gc_match/wgEncodeEH002288/DNAShape/CNN/0 | HelT motif
output/encode_101/gc_match/wgEncodeEH002288/DNAShape/CNN/1 | MGW motif
output/encode_101/gc_match/wgEncodeEH002288/DNAShape/CNN/2 | ProT motif
output/encode_101/gc_match/wgEncodeEH002288/DNAShape/CNN/3 | Roll motif

## Predict TF-DNA Binding Specicitity Using Gated-CNN (GCNN) and Long DNA Sequence
```
cd code/
python train.py --start_index 0 --end_index 1 --peak_flank 500 --network GCNN --feature_format Seq
```
```--network GCNN``` indicates that GCNN is used for model training <br/>
```--peak_flank 500``` indicates that the peak length is (2 * 500 + 1) = 1001 base pairs <br/>

### Output
The trained model and its AUC (```Test_result.txt```) on test data is located at ```output/encode_1001/gc_match/wgEncodeEH002288/Seq/GCNN```.

## Citation
If you use DESSO in your research, please cite the following paper:</br>
```
**Systematic Prediction of Regulatory Motifs from Human ChIP-Sequencing Data Based on a Deep Learning Framework,**</br>
Jinyu Yang, Adam D. Hoppe, Bingqiang Liu, Qin Ma,</br>
bioRxiv 417378; doi: https://doi.org/10.1101/417378
```
