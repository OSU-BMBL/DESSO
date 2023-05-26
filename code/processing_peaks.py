# -*- coding: utf-8 -*-
"""
Created on 08/05/2019
Updated on 10/05/2021
@author: Shuangquan Zhang
"""
import re
import os
import numpy as np
from random import shuffle
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--name',default='',help='bed')
parser.add_argument('--peak_flank',type=int, default=50, help='peaks length')

args = parser.parse_args()

def get_peaks(args):
    name=args.name
    name=name.rsplit('/',1)[1]
    name=name.split('.')[0]
    file=open(args.name,'r')
    data=file.readlines() 
    file.close()
    sortn=[]
    ll=len(data)
    for i in range(ll):
        t=data[i]
        t1=t.strip('\n').split('\t')[4]
#        print(t1)
        sortn.append(t1)
    sortn=np.array(sortn,dtype=np.float32)
    sortn_index=sortn.argsort()[::-1]
    name='../data/TfbsUniform_hg19_ENCODE/'+name+'_encode.narrowPeak'
    file1=open(name,'w')
    pattern=re.compile(r'chr\d+')
    for i in range(ll):
        t=data[sortn_index[i]]
        match_res = pattern.match(t)
        if match_res:
            t1=t.split('\t')
            tt=str(t1[0])+'\t'+ str(t1[1]) + '\t' + str(t1[2])+'\t'+str(t1[4])
            file1.writelines(tt) 
    file1.close()
    cmd='gzip '+ name
    os.system(cmd)
############################################################################

def getchrom_index(args):
    data=[]
    name=args.name
    name=name.rsplit('/',1)[1]
    name=name.split('.')[0]
    file=open(args.name,'r')
    data1=file.readlines() 
    file.close()
    sortn=[]
    ll=len(data1)
    for i in range(ll):
        t=data1[i]
        t1=t.strip('\n').split('\t')[4]
        sortn.append(t1)
    sortn=np.array(sortn,dtype=np.float32)
    sortn_index=sortn.argsort()[::-1]
    test_val=['train.bed','test.bed']
    test_val_fa=[name+'train1.fa',name+'test1.fa']
    pattern=re.compile(r'chr')
    file0=open(test_val[0],'w')
    file1=open(test_val[1],'w')
    for j in range(ll):
        t=data1[sortn_index[j]]
        match_res = pattern.match(t)
        t1=t.split('\t')
        if match_res:
            t11=round((int(t1[1])+int(t1[2]))/2)-args.peak_flank+1
            t12=round((int(t1[1])+int(t1[2]))/2)+args.peak_flank+2
            if t11 >0 and t12 > 0:
                    
                tt=str(t1[0])+'\t'+ str(int(t11)) + '\t' + str(int(t12)) +'\t'+str(t1[4])+'\n'
                if j<=1000:
                    if j % 2==0:
                        file0.writelines(tt) 
                    else:
                        file1.writelines(tt) 
                else:
                    file0.writelines(tt)
                    
    file0.close()
    file1.close()
    for k in range(2):
        cmd='bedtools getfasta -fi ../data/hg38.fa -bed '+ test_val[k]+ ' -fo '+test_val_fa[k]
        os.system(cmd)
        os.remove(test_val[k])


def makeseq(input_name,out_name):
    file=open(input_name,'r')
    data1=open(out_name,'w')
    data=file.readlines()
    file.close()
    times=0
    for t1 in data:
        if t1.startswith('>'):
            tt=0
        else:
            ttt='A'+'\t'+'peaks'+'\t'+t1.strip().upper()+'\t'+'1'+'\n'
            data1.writelines(ttt) 
            times+=1
    data1.close() 
###

def get_data(args):
    name=args.name
    name=name.rsplit('/',1)[1]
    name=name.split('.')[0]
    innames=[name+'train1.fa',name+'test1.fa']
    outnames=[name+'_encode_AC.seq',name+'_encode_B.seq',name+'_encode.seq']
    for i in range(2):
        makeseq(innames[i],outnames[i])
        cmd='cat '+outnames[i]+' >> '+outnames[2]
        os.system(cmd)
        cmd='gzip '+ outnames[i]
        os.system(cmd)
        cmd='mv '+ outnames[i]+'.gz'+' ../data/encode_'+str(args.peak_flank*2+1)
        os.system(cmd)
#        os.remove(innames[i])
    cmd='gzip '+ outnames[2]
    os.system(cmd)
    cmd='mv '+ outnames[2]+'.gz'+' ../data/encode_'+str(args.peak_flank*2+1)
    os.system(cmd)
    file=open('../data/encode_tfbs.txt','a+')
    txxt=name+'_encode'+'\t'+name+'_encode'+'\n'
    file.writelines(txxt)
    file.close()
    os.remove(innames[0])
    os.remove(innames[1])


if __name__=='__main__':
    get_peaks(args)
    getchrom_index(args)
    get_data(args)
    
            
