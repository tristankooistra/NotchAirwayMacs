import pegasus as pg
import os
import matplotlib.colors as clr
import scanpy as sc
import matplotlib.pyplot as plt

# Provide file path to h5ad data object
mnp_path = ""

# Load data
mnp_harmonized = pg.read_input(mnp_path)
mnp_harmonized = mnp_harmonized.to_anndata()

annotation_dict = {"1": "MC1 (CXCL10)",
"2": "MC2 (SPP1)",
"3": "MC3 (AREG)",
"4": "Mac1 (FABP4)",
"5": "quiesMac",
"6": "quiesMC",
"7": "Cycling (PCLAF)",
"8": "MC4 (CCR2)",
"9": "Mac2 (A2M)",
"10": "pDC (TCF4)",
"11": "migDC (CCR7)",
"12": "DC1 (CLEC9A)",
"13": "DC2 (CD1C)",
"14": "AS DC (AXL)"}

mnp_harmonized.obs["annotation"] = [annotation_dict[c] for c in mnp_harmonized.obs["new_clusters"]]


# set cluster order
cluster_order = ["12", "13", "9", "4", "14", "11", "10", "8", "1", "2", "3", "7", "6", "5"]

# associate name and cluster
cluster_name_order = [annotation_dict[x] for x in cluster_order]

# placeholder to easily reference other genes for comparison
dot_gene_dict = {
"Specialization": ["C1QC", "CSF1R", "CCR2", "CX3CR1", "ITGAX", "CIITA", "HLA-DQB1", "CD74", "MMP12", "AXL"]
}

sc.pl.violin(mnp_harmonized, keys = "ITGAX", groupby = "annotation", show = False)
plt.xlabel("")
figure = plt.gcf()
figure.set_size_inches(24, 4)

# need to enter output destination
plt.savefig("violin.pdf")

plt.show()
plt.close()