import pegasus as pg
import os
import matplotlib.colors as clr
import scanpy as sc
import matplotlib.pyplot as plt

# Change file path to h5ad data object
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

cluster_order = ["12", "13", "9", "4", "14", "11", "10", "8", "1", "2", "3", "7", "6", "5"]
cluster_name_order = [annotation_dict[x] for x in cluster_order]

# can substitute with desired gene list
dot_gene_dict = {
"NOTCH": ["DLL1", "DLL3", "JAG1", "JAG2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"]
}

dot_cmap = clr.LinearSegmentedColormap.from_list('gene_cmap', ["#d3d3d3", "#a5cde3", "#6dafd6", '#08306b'], N=200)

for key in dot_gene_dict.keys():
    plot = sc.pl.dotplot(mnp_harmonized, dot_gene_dict[key], groupby="annotation",
    categories_order=cluster_name_order,
    show=False, return_fig=True, title="{key} markers".format(key=key),
    cmap=dot_cmap, dot_max=1)
    axes_dict = plot.get_axes()
    axes_dict["mainplot_ax"].set_axisbelow(True)
    axes_dict["mainplot_ax"].grid()
    plt.subplots_adjust(left=0.2, bottom = 0.2)

    # adjust output file location
    plt.savefig("MNP_Disc_Notch_dotplot.pdf")
    plt.show()
    plt.close()