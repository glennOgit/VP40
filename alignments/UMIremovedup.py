
# cat R134A_R2.fastq | paste - - - - > R134A_R2.fastq.paste
# pipe data from bam file via samtools: samtools view -h CTD.srt.bam | python UMIgenerator.py [fastq]
# samtools view -h R134A.srt.bam | python UMIgenerator2.py R134A_R1.fastq.paste R134A_R2.fastq.paste
#Consider only sending mapped reads into the script
import pandas as pd
import numpy as np
import sys
import csv

# Takes pasted fastq file as arguments. (Need Read 1 and read 2)
fqR1df=sys.argv[1]
fqR2df=sys.argv[2]
# This loads a SAM file into a pandas dataframe. 
my_cols = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ","CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ","QUAL", "NM", "MD", "AS", "XS","SA"]
fq1_cols = ["QNAME","seq","dir","qual","tag1"]
fq2_cols = ["QNAME","seq","dir","qual","tag2"]

samdf = pd.read_csv(sys.stdin, sep='\t',header=None,names=my_cols,dtype={"MD": str,"AS": str,"XS": str,"SA": str})

hedrdf = samdf[0:0]
hedrdf = samdf[~samdf['QNAME'].astype(str).str.startswith('SN')]
hedrdf = hedrdf[hedrdf.QNAME != '@PG']
hedrdf = hedrdf.dropna(axis=1,how='all')

samdf = samdf[samdf['QNAME'].astype(str).str.startswith('SN')]
samdf = samdf[samdf.TLEN != 0] #deleting unaligned reads where TLEN equal 0
samdf = samdf.sort_values(by=['QNAME'],axis=0)
samdf[['PNEXT', 'TLEN']] = samdf[['PNEXT', 'TLEN']].astype(int)

fqR1df = pd.read_csv(fqR1df, sep='\t',header=None,names=fq1_cols) # get read 1 file info
fqR1df['QNAME'] = fqR1df['QNAME'].str[1:-13] # this drops first char of string and last 12 chars of string
fqR1df["tag1"] = fqR1df["seq"].astype(str).str[:4] # gets first 4nt of read and sets creates its own column
fqR1df = fqR1df[['QNAME','tag1']]

fqR2df = pd.read_csv(fqR2df, sep='\t',header=None,names=fq2_cols) # get read 2 file info
fqR2df['QNAME'] = fqR2df['QNAME'].str[1:-13] # This edits qname to match bam file
fqR2df["tag2"] = fqR2df["seq"].astype(str).str[:4] # gets first 4nt of read and sets creates its own column
fqR2df = fqR2df[['QNAME','tag2']]


mergedf = pd.merge(samdf, fqR1df, how="left", on="QNAME")
finaldf = pd.merge(mergedf, fqR2df, how="left", on="QNAME")
finaldf['tag'] = finaldf['tag1'].astype(str) + finaldf['tag2']
finaldf = finaldf.drop(['tag1','tag2'], axis=1)
finaldf = finaldf.sort_values(by=['tag','RNAME','POS'],axis=0)


tempdf = finaldf.copy()
tempdf = tempdf.drop_duplicates(subset=["QNAME","RNAME"])
tempdf = tempdf.sort_values(by=['tag','MAPQ','RNAME','POS'],axis=0,ascending=False)

tempdf['dups'] = tempdf.duplicated(subset=["tag","RNAME","POS"])
tempdf = tempdf[tempdf.dups]

toremove = tempdf["QNAME"].tolist()
finaldf = finaldf[~finaldf['QNAME'].isin(toremove)]
finaldf = finaldf.drop('tag',axis=1)
finaldf[['PNEXT', 'TLEN']] = finaldf[['PNEXT', 'TLEN']].astype(int)

hedrdf.to_csv(sys.stdout,sep='\t',header=False,index=False)
finaldf.to_csv(sys.stdout,sep='\t',header=False,index=False)


