#code for comparing the DupGenFinder, doubletrouble, DupyliCate results on A. thaliana Col-0 dataset

# final benchmarking code without BUSCO results of DupyliCate

import pandas as pd

# Read the TSV file
dft = pd.read_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_tandems.tsv", sep="\t")

# Extract the first column, split by commas, flatten into one list
genest = []
for entry in dft.iloc[:, 0].dropna():
    for gene in entry.split(","):
        genest.append(gene.strip())
        
# Extract the first column (excluding header) and save
pd.Series(genest).to_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_tandems_processed.tsv", index=False, header=False, sep="\t")

# Read the TSV file
dfp = pd.read_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_proximals.tsv", sep="\t")

# Extract the first column, split by commas, flatten into one list
genesp = []
for entry in dfp.iloc[:, 0].dropna():
    for gene in entry.split(","):
        genesp.append(gene.strip())
        
# Extract the first column (excluding header) and save
pd.Series(genesp).to_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_proximals_processed.tsv", index=False, header=False, sep="\t")

# Read the TSV file
dfd = pd.read_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_dispersed_duplicates.tsv", sep="\t")

# Extract the first column, split by commas, flatten into one list
genesd = []
for entry in dfd.iloc[:, 0].dropna():
    for gene in entry.split(","):
        genesd.append(gene.strip())
        
# Extract the first column (excluding header) and save
pd.Series(genesd).to_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_dispersed_duplicates_processed.tsv", index=False, header=False, sep="\t")

# Read the TSV file
dfs = pd.read_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_singletons.tsv", sep="\t")

# Extract the first column, split by commas, flatten into one list
geness = []
for entry in dfs.iloc[:, 0].dropna():
    for gene in entry.split(","):
        geness.append(gene.strip())
        
# Extract the first column (excluding header) and save
pd.Series(geness).to_csv("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_singletons_processed.tsv", index=False, header=False, sep="\t")


with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/Dupgenfinder_Results/dupgenfinder_proximals.tsv",'r') as f1:
    dupgenfinder_tandems_list = f1.read().splitlines()
    dupgenfinder_tandems_set = set (dupgenfinder_tandems_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/Dupgenfinder_Results/dupgenfinder_tandems.tsv",'r') as f2:
    dupgenfinder_proximals_list = f2.read().splitlines()
    dupgenfinder_proximals_set = set (dupgenfinder_proximals_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/Dupgenfinder_Results/dupgenfinder_misc_dups.tsv",'r') as f3:
    dupgenfinder_misc_list = f3.read().splitlines()
    dupgenfinder_misc_set = set (dupgenfinder_misc_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/Dupgenfinder_Results/dupgenfinder_singletons.tsv",'r') as f4:
    dupgenfinder_singletons_list = f4.read().splitlines()
    dupgenfinder_singletons_set = set (dupgenfinder_singletons_list)
print('no. of dupgenfinder tandem genes are '+str(len(dupgenfinder_tandems_set)))
print('no. of dupgenfinder proximal genes are '+str(len(dupgenfinder_proximals_set)))
print('no. of dupgenfinder misc genes are '+str(len(dupgenfinder_misc_set)))
print('no. of dupgenfinder singletons genes are '+str(len(dupgenfinder_singletons_set)))
print ('\n')

with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/doubletrouble_results/doubletrouble_proximals.tsv",'r') as f1:
    doubletrouble_tandems_list = f1.read().splitlines()
    doubletrouble_tandems_set = set (doubletrouble_tandems_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/doubletrouble_results/doubletrouble_tandems.tsv",'r') as f2:
    doubletrouble_proximals_list = f2.read().splitlines()
    doubletrouble_proximals_set = set (doubletrouble_proximals_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/doubletrouble_results/doubletrouble_misc_dups.tsv",'r') as f3:
    doubletrouble_misc_list = f3.read().splitlines()
    doubletrouble_misc_set = set (doubletrouble_misc_list)
print('no. of doubletrouble tandem genes are '+str(len(doubletrouble_tandems_set)))
print('no. of doubletrouble proximal genes are '+str(len(doubletrouble_proximals_set)))
print('no. of doubletrouble misc genes are '+str(len(doubletrouble_misc_set)))
print ('\n')

with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_tandems_processed.tsv",'r') as f1:
    dupylicate_tandems_list = f1.read().splitlines()
    dupylicate_tandems_set = set (dupylicate_tandems_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_proximals_processed.tsv",'r') as f2:
    dupylicate_proximals_list = f2.read().splitlines()
    dupylicate_proximals_set = set (dupylicate_proximals_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_dispersed_duplicates_processed.tsv",'r') as f3:
    dupylicate_disp_list = f3.read().splitlines()
    dupylicate_disp_set = set (dupylicate_disp_list)
with open("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/dupylicate_v1.0_results_for_benchmarking_final/Without_busco/AthCol0_singletons_processed.tsv",'r') as f4:
    dupylicate_singletons_list = f4.read().splitlines()
    dupylicate_singletons_set = set (dupylicate_singletons_list)
print('no. of dupylicate tandem genes are '+str(len(dupylicate_tandems_set)))
print('no. of dupylicate proximal genes are '+str(len(dupylicate_proximals_set)))
print('no. of dupylicate dispersed genes are '+str(len(dupylicate_disp_set)))
print('no. of dupylicate singletons genes are '+str(len(dupylicate_singletons_set)))
print ('\n')

tandprototset_dupgenfinder = dupgenfinder_tandems_set | dupgenfinder_proximals_set
print('tand+prox dupgenfinder ' + str(len(tandprototset_dupgenfinder)))

tandprototset_doubletrouble = doubletrouble_tandems_set | doubletrouble_proximals_set
print('tand+prox doubletrouble ' + str(len(tandprototset_doubletrouble)))

tandprototset_dupylicate = dupylicate_tandems_set | dupylicate_proximals_set
print('tand+prox dupylicate ' + str(len(tandprototset_dupylicate)))
print ('\n')


dupstotset_dupgenfinder = dupgenfinder_tandems_set | dupgenfinder_proximals_set | dupgenfinder_misc_set
print('tot dups dupgenfinder ' + str(len(dupstotset_dupgenfinder)))

dupstotset_doubletrouble = doubletrouble_tandems_set | doubletrouble_proximals_set | doubletrouble_misc_set
print('tot dups doubletrouble ' + str(len(dupstotset_doubletrouble)))

dupstotset_dupylicate = dupylicate_tandems_set | dupylicate_proximals_set | dupylicate_disp_set
print('tot dups dupylicate ' + str(len(dupstotset_dupylicate)))

print ('\n')

common = dupstotset_dupgenfinder & dupstotset_doubletrouble & dupstotset_dupylicate
print(len(common))
print('\n')
i=0
#case a
elementsa=[]
for ele in dupstotset_dupylicate:
    if ele not in dupstotset_doubletrouble | dupstotset_dupgenfinder:
        i+=1
        elementsa.append(ele)
print('duplicates in dupylicate but not shared with doubletrouble and dupgenfinder '+str(i)+' in number and are'+ str(elementsa))
j=0
#case b
elementsb=[]
for ele in dupstotset_doubletrouble | dupstotset_dupgenfinder:
    if ele not in dupstotset_dupylicate:
        j+=1
        elementsb.append(ele)
print('duplicates in doubletrouble and dupgenfinder but not shared with dupylicate '+str(j)+' in number and are'+ str(elementsb))
k=0
#case c
elementsc=[]
for ele in tandprototset_dupylicate:
    if ele not in tandprototset_doubletrouble | tandprototset_dupgenfinder:
        k+=1
        elementsc.append(ele)
print('tand prox duplicates in dupylicate but not shared with doubletrouble and dupgenfinder tand prox duplicates '+str(k)+' in number and are'+ str(elementsc))
l=0
#case d
elementsd=[]
for ele in tandprototset_doubletrouble | tandprototset_dupgenfinder:
    if ele not in tandprototset_dupylicate:
        l+=1
        elementsd.append(ele)
print('tand prox duplicates in doubletrouble and dupgenfinder but not shared with dupylicate tand prox duplicates '+str(l)+' in number and are'+ str(elementsd))

# NOTE: The same code above was used with the results of DupyliCate’s BUSCO-based thresholding method to do the benchmarking

# Runtime benchmarking of DupyliCate, DupGen_finder and doubletrouble

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Data
data = {
    "Tools": [
        "DupyliCate (with BUSCO)",
        "DupyliCate (without BUSCO)",
        "doubletrouble",
        "DupGen_finder"
    ],
    "Run1": [14.39438, 9.14326, 6.78565, 13.8952],
    "Run2": [14.223226, 9.11478, 6.75687, 13.88117],
    "Run3": [14.2873456, 9.12103, 6.82438, 13.88183],
    "Run4": [14.2169789, 9.147861, 6.74665, 13.8717],
    "Run5": [14.2517465, 9.135745, 6.756, 13.871766]
}

df = pd.DataFrame(data)
df.set_index("Tools", inplace=True)

# Mean and std
means = df.mean(axis=1)
stds = df.std(axis=1)

plt.figure(figsize=(9, 6))
x = np.arange(len(means))

# Plot thinner bars with error bars
bars = plt.bar(
    x, means.values, yerr=stds.values,
    capsize=5, width=0.4, color="lightgreen", edgecolor="black", label="Mean ± SD"
)

# Add replicate points (jittered for visibility)
for i, tool in enumerate(df.index):
    y = df.loc[tool].values
    jitter = np.random.uniform(-0.1, 0.1, size=len(y))
    plt.scatter(np.full_like(y, i) + jitter, y,
                color="gray", alpha=0.7, s=40, label="Replicate runs" if i == 0 else "")

# Annotate n above each bar
for i, bar in enumerate(bars):
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, height + 0.2, "n=5",
             ha='center', va='bottom', fontsize=10)

plt.xticks(x, means.index, rotation=20, ha="right")
plt.ylabel("Runtime (minutes)")

# Adjust y-axis to create space at the top
plt.ylim(0, 18)

# Better legend positioning with transparency
plt.legend(loc="upper right", framealpha=0.95, bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()
plt.savefig("/home/shakunthala/Downloads/mpsBigData/Shakunthala_Uni_Bonn_Backup/PhD_Project/Gene_Duplications/Dupylicate/Final_dupylicate_figures/V2/Benchmark_runtimes.png",dpi=600)
plt.show()