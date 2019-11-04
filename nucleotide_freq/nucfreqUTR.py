
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pybedtools
from pybedtools import BedTool
import seaborn as sns; sns.set()

#Figure out how to get strand orientation of each gene. 


df = pd.read_csv('3UTR.sense2.u.fa', sep='\t',names='s')
df = df[~df['s'].astype(str).str.startswith('>')]
df['s'] = df['s'].str.upper()

df = df['s'].apply(lambda x: pd.Series(list(x)))
df = df.iloc[:,0:400]
print(df)
df2 = pd.DataFrame({'A' : []})
for column in df.columns:
    df2 = df2.append(df[column].value_counts())
df2 = df2.divide(9588)

print(df2)



sns.set_style("darkgrid")
plt.figure(figsize=(10,3))
lm = plt.plot(df2['A'],color='limegreen')
plt.plot(df2['T'],color='red')
plt.plot(df2['C'],color='blue')
plt.plot(df2['G'],color='gold')
plt.ylim( (0.08, .6) )

plt.legend(['A', 'U', 'C', 'G'], loc='upper left')
plt.savefig('3UTR.sense2.u.png',dpi=300,figsize=(20,6))