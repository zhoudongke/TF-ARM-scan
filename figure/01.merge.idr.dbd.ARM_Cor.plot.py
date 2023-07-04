#!/usr/bin/env python
# coding: utf-8

# In[1]:

import argparse
from scipy.signal import correlate
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import math
from matplotlib.ticker import FormatStrFormatter
import os


# In[2]:
## parse the arguments ##

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--family', type=str, help='family index file')
parser.add_argument('-i', '--idr', type=str, help='IDRs score file(eg. file.csv)')
parser.add_argument('-d', '--dbd', type=str, help='DBD region file, with 6 columns')
parser.add_argument('-a', '--fasta', type=str, help='TFs fasta file')

args = parser.parse_args()

family_index = args.family
idr_file = args.idr
dbd_file = args.dbd
tf_fasta_file = args.fasta


def digitized(seq,start,tf):
    digitized_seq = ""
    for i,aa in enumerate (seq):
        if i+start in idr_index[tf] and float(idr_index[tf][i+start]) < 0.2:
            digitized_seq += "0"
        elif aa in ["R", "K"]:
            digitized_seq += "1"
        else:
            digitized_seq += "0"
    if 'R' not in seq:
        digitized_seq = "000000000"
    return digitized_seq


# In[3]:


family={}
with open (family_index) as file:
    for line in file:
        line=re.sub("\n$","",line)
        tf=re.split("\t",line)[0]
        fa=re.split("\t",line)[1]
        if fa not in family:
            family[fa]=[]
            family[fa].append(tf)
        else:
            family[fa].append(tf)


# In[4]:


dbd={}
dbd_df=pd.read_table(dbd_file,
                     header=None,
                     index_col=0,
                     usecols=[0, 2, 3]
                    )
#dbd_df
dbd=dbd_df.to_dict('index')  #只有to_dict()的第一个参数用“index”，才会生成以行名作为索引的字典。要以列明作为索引，要用“columns”为第一个参数（orient）


# In[5]:


idr_index={}
with open(idr_file) as idr:
    for lines in idr:
        tf_name=re.split("\s",lines.split(',')[1])[1]
        for i,aa in enumerate(lines.split(',')[3:]):
            if tf_name in idr_index:
                idr_index[tf_name][i]=aa
            else:
                idr_index[tf_name]={}
                idr_index[tf_name][i]=aa


# In[6]:


'''
search_kernel = "111110111"
digitized_seq = "111110111"
correlation = correlate(
                    [int(d) for d in search_kernel],
                    [int(d) for d in digitized_seq],
                    method="direct",
                    mode='valid'
                    )
print (correlation)
'''


# In[32]:


search_kernel = "111110111"
plot={}
with open(tf_fasta_file) as file:
    content = file.read()
    for block in content.split(">"):
        if block:
            lines = block.split("\n")
            lines[0] = re.sub("^>","",lines[0])
            tf_name = re.split("\s+",lines[0])[0]
            sequence = "".join(lines[1:])
            length=len(sequence)
            #print(length)
            if tf_name not in idr_index or tf_name not in dbd:
                continue
            for i in range(length-8):
                digitized_seq = digitized(seq=sequence[i:i+9],start=i,tf=tf_name)
                correlation = correlate(
                    [int(d) for d in search_kernel],
                    [int(d) for d in digitized_seq],
                    method="direct",
                    mode='valid'
                    )
                pos=i+4
                if tf_name not in plot:
                    #plot[tf_name] = {'Pos':[], 'DBD':[], 'Cor':[]}
                    plot[tf_name] = []
                #plot[tf_name]['Pos'].append(pos)
                if pos >= dbd[tf_name][2] and pos <= dbd[tf_name][3]:
                    plot[tf_name].append({'Pos':pos , 'DBD':0.2 , 'Cor':correlation[0]/8, 'IDR' : float(idr_index[tf_name][pos])})
                else:
                    plot[tf_name].append({'Pos':pos , 'DBD':0 , 'Cor':correlation[0]/8 , 'IDR' : float(idr_index[tf_name][pos])})


# In[36]:


'''
fig,ax=plt.subplots()
df = pd.DataFrame(plot['AT5G17810.2'])
ax.set_ylim([0, 1])
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.bar(df['Pos'],df['Cor'],color='#000000')
ax.bar(df['Pos'],df['DBD'],color='#cc9999')
ax.plot(df['Pos'],df['IDR'],color='#6699cc',alpha=0.5)
ax.set_title(tf)
'''


# In[37]:


'''
new_d = {k: v for k, v in plot.items() if k in family['WOX']}
fig,axs=plt.subplots(math.ceil(len(new_d.keys())/2),2,
                    figsize=(12,len(new_d.keys()))
                    )
for tf,ax in zip(new_d.keys(),axs.flatten()):
#print (tf)
    df = pd.DataFrame(plot[tf])
    ax.set_ylim([0, 1])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.bar(df['Pos'],df['Cor'],color='#000000')
    ax.bar(df['Pos'],df['DBD'],color='#cc9999')
    ax.plot(df['Pos'],df['IDR'],color='#6699cc',alpha=0.5)
    ax.set_title(tf)
plt.tight_layout()
'''


# In[10]:

if not os.path.exists("family.plot"):
    os.makedirs("family.plot")


for fa in family.keys():
    new_d = {k: v for k, v in plot.items() if k in family[fa]}
    if len(new_d.keys()) ==0:
        continue
    fig,axs=plt.subplots(math.ceil(len(new_d.keys())/2),2,
                        figsize=(6,len(new_d.keys())/2)
                        )
    for tf,ax in zip(new_d.keys(),axs.flatten()):
    #print (tf)
        df = pd.DataFrame(plot[tf])
        ax.set_ylim([0, 1])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.bar(df['Pos'],df['Cor'],color='#000000')
        ax.bar(df['Pos'],df['DBD'],color='#cc9999')
        ax.plot(df['Pos'],df['IDR'],color='#6699cc',alpha=0.5)
        ax.set_title(tf)
    plt.tight_layout()
    fa=re.sub('\/','_',fa)
    plt.savefig(f"family.plot/{fa}.pdf")
    plt.close(fig)
    plt.ioff()

