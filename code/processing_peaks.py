# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 21:05:21 2019

@author: Shuangquan Zhang, Cankun Wang
"""
import re
import os
import numpy as np
from random import shuffle
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--name',default='',help='bed')
args = parser.parse_args()

def get_peaks(args):
	#strname=name.split('_')[0]
	#name="../data/Fox01/fox01_peaks.bed.gz"
	name=args.name
	name=name.split('.')[2].split('/')[3]
	file=open(args.name,'r')
	data=file.readlines() 
	file.close()
	name='../data/TfbsUniform_hg19_ENCODE/'+name+'_encode.narrowPeak'
	file1=open(name,'w')
	pattern=re.compile(r'chr\d+')
	for t in data:
		match_res = pattern.match(t)
		if match_res:
			t1=t.split('\t')
			tt=str(t1[0])+'\t'+ str(t1[1]) + '\t' + str(t1[2])+'\t'+str(t1[4])
#                tt=tt.replace('t','    ')
			file1.writelines(tt) 
	file1.close()
	cmd='gzip '+ name
	os.system(cmd)


def sh_str(s):
    str_l=list(s.strip('\n'))
    shuffle(str_l)
    strs=''
    return strs.join(str_l)

def getchrom_index(args):
    
    file=open(args.name,'r')
    data=file.readlines() 
    ll=len(data)
    file.close()
    test_val=['train.bed','test.bed']
    test_val_fa=['train1.fa','test1.fa']
    train_index=round(ll*0.8)
    test_val_index=np.array([[0,train_index],[train_index,ll]],dtype=np.int32)
    pattern=re.compile(r'chr\d+')
    for i in range(2):
        file1=open(test_val[i],'w')
        for j in range(test_val_index[i][0],test_val_index[i][1]):
            t=data[j]
            match_res = pattern.match(t)
            if match_res:
                t1=t.split('\t')
                t11=round((int(t1[1])+int(t1[2]))/2)-499
                t12=round((int(t1[1])+int(t1[2]))/2)+502
                tt=str(t1[0])+'\t'+ str(int(t11)) + '\t' + str(int(t12)) +'\n'
                file1.writelines(tt) 
        file1.close()
        cmd='bedtools getfasta -fi ../data/hg38.fa -bed '+ test_val[i]+ ' -s -fo '+test_val_fa[i]
        os.system(cmd)
        os.remove(test_val[i])


def makeseq(input_name,out_name):
    file=open(input_name,'r')
    data1=open(out_name,'w')
    data=file.readlines()
    file.close()
    for t1 in data:
        if t1.startswith('>'):
            tt=0
        else:
            ttt='A'+'\t'+'peaks'+'\t'+t1.strip().upper()+'\t'+'1'+'\n'
            data1.writelines(ttt) 
            ttt=''
    data1.close() 


def get_data(args):
    name=args.name
    name=name.split('.')[2].split('/')[3]
    innames=['train1.fa','test1.fa']
    outnames=[name+'_encode_AC.seq',name+'_encode_B.seq',name+'_encode.seq']
    for i in range(2):
        makeseq(innames[i],outnames[i])
        cmd='cat '+outnames[i]+' >> '+outnames[2]
        os.system(cmd)
        cmd='gzip '+ outnames[i]
        os.system(cmd)
        cmd='mv '+ outnames[i]+'.gz'+' ../data/encode_1001'
        os.system(cmd)
        os.remove(innames[i])
    cmd='gzip '+ outnames[2]
    os.system(cmd)
    cmd='mv '+ outnames[2]+'.gz'+' ../data/encode_1001'
    os.system(cmd)
    with open('../encode_tfbs.txt', 'a') as the_file:
        the_file.write(name + '_encode\t' + name + '_encode')
        
        
        
'''def w_tfbs(args):
    file=open('../data/encode_tfbs.txt','w')
    name=args.name
    name=name.split('.')[0]+'_wgEncodeAwgTfbsBroadDnd41CtcfUniPk'+'\t'+name.split('.')[0]+'_wgEncodeAwgTfbsBroadDnd41CtcfUniPk'
    fiile.writelines(name)
    '''
    
    
    


if __name__=='__main__':
#    makeseq('peaks1001.fa')
    get_peaks(args)
    getchrom_index(args)
    get_data(args)
    
            
