import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

args = sys.argv

family=pd.read_table(args[1],header=None)

plt.figure(figsize=(15,5))
ax=sns.barplot(family,x=0,y=1,color="#0099cc")
ax.set_xlabel("TF family", fontsize=14)
ax.set_ylabel("counts", fontsize=14)
plt.xticks(
    rotation=45, 
    horizontalalignment='right',
    #fontweight='light',
    #fontsize='x-large'  
)
plt.savefig(f"family.plot/TF.family.counts.barplot.pdf")
plt.savefig(f"family.plot/TF.family.counts.barplot.png")

plt.figure(figsize=(15,5))
ax=sns.barplot(family,x=0,y=2,color="#0099cc")
ax.set_xlabel("TF family", fontsize=14)
ax.set_ylabel("Ratio_of_family", fontsize=14)
plt.axhline(y=0.54, color='r', linestyle='--')
plt.xticks(
    rotation=45, 
    horizontalalignment='right',
    #fontweight='light',
    #fontsize='x-large'  
)
plt.savefig(f"family.plot/TF.family.ratio.barplot.pdf")
plt.savefig(f"family.plot/TF.family.ratio.barplot.png")