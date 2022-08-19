import os
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors, rc, rcParams
from numpy import log
from scipy.stats import pearsonr

import greatpy as great

plt.rcParams.update({"font.size": 14, "font.weight": "normal"})


def scatterplot(
    great_df: pd.DataFrame,
    colname_x: str,
    colname_y: str,
    title: str = "",
    minus_log10=True,
    ax: matplotlib.axes or None = None,
) -> None:
    """
    Create a scatterplot from a
    pandas dataframe between two columns.\n
    A logarithmic scale can be used.\n

    Parameters
    ----------
    great_df : pd.DataFrame
        Output of the greatpy.tl.enrichment function\n
    colname_x : str
        Name of the column to be used as x axis\n
    colname_y : str
        Name of the column to be used as y axis\n
    title : str
        Title of the plot.
        Default is `""`\n
    minus_log10 : bool
        If True, the logarithmic scale is used\n
        Default is `True`\n
    ax : matplotlib.axes or None
        Define the position of the plot in a figure\n
        Default is `None`\n

    Returns
    -------
    None
        Scatterplot between the two columns passed in args

    """
    great_df = great_df.dropna()
    great_df = great_df.astype({colname_x: "float64", colname_y: "float64"})
    if minus_log10:
        great_df[f"-log({colname_x})"] = -log(great_df[colname_x])
        great_df[f"-log({colname_y})"] = -log(great_df[colname_y])
        g = sns.scatterplot(data=great_df, x=f"-log({colname_x})", y=f"-log({colname_y})", ax=ax)
        g.set_title(title, fontsize=18)
        g.set_xlabel(f"-log({colname_x})", fontsize=12)
        g.set_ylabel(f"-log({colname_y})", fontsize=12)
    else:
        g = sns.scatterplot(data=great_df, x=colname_x, y=colname_y, ax=ax)
        g.set_title(title, fontsize=20)
        g.set_xlabel(colname_x, fontsize=15)
        g.set_ylabel(colname_y, fontsize=15)


def graph_nb_asso_per_peaks(
    test: str or pd.DataFrame,
    regdom: str or pd.DataFrame,
    ax: matplotlib.axes or None = None,
    color: str = None,
) -> None:
    """
    Creates a barplot representing the
    percentage of peaks for all possible association numbers

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested\n
    regdom : str or pd.DataFrame
        Regulatory domain of all genes in the genome\n
    ax : matplotlib.axes or None
        Define the position of the plot in a figure\n
        Default is `None`\n
    color : str
        Color of the bar\n
        Default is `None` \n

    Returns
    -------
    None
        Barplot of the number of associated genes per peak

    """
    nb_asso_per_peaks = great.tl.get_nb_asso_per_region(test, regdom)

    nb = {
        "number": [],
        "number_genes": [],
        "percentage": [],
    }
    for i in list(set(nb_asso_per_peaks.values())):
        nb["number"].append(i)
        nb["number_genes"].append(list(nb_asso_per_peaks.values()).count(i))
        nb["percentage"].append(
            round((list(nb_asso_per_peaks.values()).count(i) / len(nb_asso_per_peaks.keys())) * 100)
        )
    nb = pd.DataFrame(nb, columns=["number", "number_genes", "percentage"], index=nb["number"])

    g = sns.barplot(data=nb, x="number", y="percentage", ax=ax, color=color)
    g.set_title("Number of associated genes per region", fontsize=20)
    g.set_xlabel("Number of associated genes per region", fontsize=15)
    g.set_ylabel("Genomic region (%)", fontsize=15)

    for i in range(nb.shape[0]):
        x = nb.iloc[i]["number"]
        y = nb.iloc[i]["percentage"]
        g.text(x=x - 0.06, y=y + 1, s=nb.number_genes[i])


def graph_dist_tss(
    test: str or pd.DataFrame,
    regdom: str or pd.DataFrame,
    ax: matplotlib.axes or None = None,
    color: str = "#325fa8",
) -> None:
    """
    Creation of a barplot of the distance
    between the peaks and the TSS of the associated gene(s).\n

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested\n
    regdom : str or pd.DataFrame
        Regulatory domain of all genes in the genome\n
    ax : matplotlib.axes or None
        Define the position of the plot in a figure\n
        Default is `None`\n
    color : str
        Color of the bar\n
        Default is `#325fa8`\n

    Returns
    -------
    None
        Barplot of the distance between the peaks and the TSS of the associated gene(s)

    """
    res = {
        "<-500": [0],
        "-500:-50": [0],
        "-50:-5": [0],
        "-5:0": [0],
        "0:5": [0],
        "5:50": [0],
        "50:500": [0],
        ">500": [0],
    }
    nb = 0

    dist = great.tl.get_dist_to_tss(test, regdom)
    for i in dist.values():
        for j in i:
            if j < -500000:
                res["<-500"][0] += 1
            elif j < -50000:
                res["-500:-50"][0] += 1
            elif j < -5000:
                res["-50:-5"][0] += 1
            elif j < 0:
                res["-5:0"][0] += 1
            elif j < 5000:
                res["0:5"][0] += 1
            elif j < 50000:
                res["5:50"][0] += 1
            elif j < 500000:
                res["50:500"][0] += 1
            else:
                res[">500"][0] += 1
            nb += 1
    df = pd.DataFrame(res).transpose().rename(columns={0: "count"})
    df["percentage"] = (df["count"] / nb) * 100
    df = df.reset_index(drop=False).rename(columns={"index": "distance"})
    g = sns.barplot(data=df, x="distance", y="percentage", color=color, ax=ax)
    for idx, p in enumerate(g.patches):
        g.annotate(str(df.iloc[idx]["count"]), (p.get_x() + p.get_width() / 2, p.get_height()))
    g.set_xlabel("Distance to TSS (kb)", fontsize=15)
    g.set_ylabel("Genomic region (%)", fontsize=15)
    g.set_title("Binned by absolute distance to TSS", fontsize=20)


def graph_absolute_dist_tss(
    test: str or pd.DataFrame,
    regdom: str or pd.DataFrame,
    ax: matplotlib.axes or None = None,
    color: str = "#325fa8",
) -> None:
    """
    Creation of a barplot of the absolute
    distance between the peaks and the TSS of the associated gene(s).

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested\n
    regdom : str or pd.DataFrame
        Regulatory domain of all genes in the genome\n
    ax : matplotlib.axes or None
        Define the position of the plot in a figure\n
        Default is `None`\n
    color : str
        Color of the bar\n
        Default is `#325fa8`\n

    Returns
    -------
    None
        Barplot of the absolute distance between the peaks and the TSS of the associated gene(s)

    """
    res = {"0:5": [0], "5:50": [0], "50:500": [0], ">500": [0]}
    nb = 0
    dist = great.tl.get_dist_to_tss(test, regdom)

    for i in dist.values():
        for j in i:
            j = abs(j)
            if j < 5000:
                res["0:5"][0] += 1
            elif j < 50000:
                res["5:50"][0] += 1
            elif j < 500000:
                res["50:500"][0] += 1
            else:
                res[">500"][0] += 1
            nb += 1
    df = pd.DataFrame(res).transpose().rename(columns={0: "count"})
    df["percentage"] = (df["count"] / nb) * 100
    df = df.reset_index(drop=False).rename(columns={"index": "distance"})
    g = sns.barplot(data=df, x="distance", y="percentage", color=color, ax=ax)
    for idx, p in enumerate(g.patches):
        g.annotate(str(df.iloc[idx]["count"]), (p.get_x() + p.get_width() / 2, p.get_height()))
    g.set_xlabel("Absolute distance to TSS (kb)", fontsize=15)
    g.set_ylabel("Genomic region (%)", fontsize=15)
    g.set_title("Binned by absolute distance to TSS", fontsize=20)


def get_all_comparison(
    results: dict,
    out_dir: str = "../data/tests/test_data/output/",
    information_folder: str = "../data/human/",
    good_gene_associations: bool = True,
    disp_scatterplot: bool = True,
    stats: bool = True,
):
    """
    Plot the comparaison between greatpy and Great from some files compute by `great.tl.enrichment_multiple`.

    Parameters
    ----------
    results : dict
        Dictionary of results from `great.tl.enrichment_multiple`\n
    out_dir : str
        Path of the output directory with the results of great webserver. \n
        Default is `../data/tests/test_data/output/`\n
    information_folder : str
        path of the folder with the information files for the tests. \n
        Default is `../data/human/`\n
        The input folder should contains the files : \n
        - `information_folder/assembly_eg_hg38/regulatory_domain.bed` \n
        - `information_folder/assembly_eg_hg38/chr_size.bed`\n
    good_gene_associations : bool
        If True, the function return the number of good gene associations\n
    disp_scatterplot : bool
        If True, the function display the scatterplot of the comparaison\n
    stats : bool
        If True, the function return the statistics of the comparaison\n

    Returns
    -------
    pp : pd.DataFrame
        Dataframe of the number of row lost between before preprocessing and after preprocessing\n
    asso : pd.DataFrame
        DataFrame of the number of good gene associations for each file\n
    stats_df : pd.DataFrame
        DataFrame of the statistics of the comparaison for each file\n

    Example
    -------
    >>> test = [
    ...    '../data/tests/test_data/input/09_ERF.bed', '../data/tests/test_data/input/10_MAX.bed',
    ...    '../data/tests/test_data/input/01_random.bed', '../data/tests/test_data/input/04_ultra_hg38.bed',
    ...    '../data/tests/test_data/input/02_srf_hg38.bed', '../data/tests/test_data/input/08_FOXO3.bed',
    ...    '../data/tests/test_data/input/06_height_snps_hg38.bed'
    ...     ]
    >>> results = great.tl.enrichment_multiple(
    ...    tests = t,
    ...    regdom_file=regdom,
    ...    chr_size_file=size,
    ...    annotation_file="../data/human/ontologies.csv",
    ...    annpath=None,
    ...    binom=True,
    ...    hypergeom=True
    ...    )
    >>> pp,asso,stat = get_all_comparison(results)

    >>> pp
    ...    |    | name   |   before_pp_greatpy_size |   before_pp_great_size |   final_size |   %_of_GO_from_great_lost |
    ...    |---:|:-------|-------------------------:|-----------------------:|-------------:|--------------------------:|
    ...    |  0 | ERF    |                     6014 |                   2410 |         1833 |                     23.94 |
    ...    |  1 | MAX    |                     2996 |                   2395 |         1481 |                     38.16 |
    ...    |  2 | random |                      579 |                    197 |          117 |                     40.61 |
    ...    |  3 | ultra  |                     3265 |                   2175 |         1393 |                     35.95 |
    ...    |  4 | srf    |                     4810 |                   2681 |         1854 |                     30.85 |

    >>> asso
    ...    |    | name   |   number_good_gene_asso |   number_genes_asso_lost |   number_gene_asso_excess |
    ...    |---:|:-------|------------------------:|-------------------------:|--------------------------:|
    ...    |  0 | ERF    |                    1456 |                        0 |                        36 |
    ...    |  1 | MAX    |                     428 |                        0 |                         4 |
    ...    |  2 | random |                      57 |                        0 |                         1 |
    ...    |  3 | ultra  |                     496 |                        0 |                         2 |
    ...    |  4 | srf    |                     923 |                        0 |                         7 |

    >>> stat
    ...    |    | name   |   pearson_binom |   pearson_hypergeom |
    ...    |---:|:-------|----------------:|--------------------:|
    ...    |  0 | ERF    |        0.57644  |            0.609552 |
    ...    |  1 | MAX    |        0.601492 |            0.670499 |
    ...    |  2 | random |        0.240765 |            0.124707 |
    ...    |  3 | ultra  |        0.52949  |            0.675438 |
    ...    |  4 | srf    |        0.631909 |            0.597787 |
    """
    pp = {
        "name": [],
        "before_pp_greatpy_size": [],
        "before_pp_great_size": [],
        "final_size": [],
        "%_of_GO_from_great_lost": [],
    }
    asso = {
        "name": [],
        "number_good_gene_asso": [],
        "number_genes_asso_lost": [],
        "number_gene_asso_excess": [],
    }

    stat_df = {"name": [], "pearson_binom": [], "pearson_hypergeom": []}

    for path in results.keys():
        id = path.split("/")[-1].split("_")[0]
        name = path.split("/")[-1].split("_")[1].split(".")[0]
        i = 0
        great_out = ""
        great_asso = ""

        for out_path in os.listdir(out_dir):
            if out_path.split("_")[0] == id:
                if re.match(".*hg19.*", out_path) is not None:
                    assembly = "hg19"
                else:
                    assembly = "hg38"

                if re.match(".*output.*", out_path) is not None:
                    great_out = out_dir + out_path
                else:
                    great_asso = out_dir + out_path

        test = path
        regdom = f"{information_folder}{assembly}/regulatory_domain.bed"

        if great_out == "" or great_asso == "":
            return False

        pp["name"].append(name)
        enrichment_tot = results[test]
        enrichment_tot = great.tl.set_bonferroni(enrichment_tot, 0.05)
        enrichment_tot = great.tl.set_fdr(enrichment_tot, 0.05)

        great_webserver = pd.read_csv(
            great_out,
            sep="\t",
            comment="#",
            names=[
                "ontologie",
                "term_name",
                "ID",
                "binom_p_value",
                "binom_bonferroni",
                "binom_fdr",
                "hyper_p_value",
                "hyper_bonferroni",
                "hyper_fdr",
            ],
            index_col=False,
            dtype={
                "term_name": "object",
                "ID": "object",
                "binom_p_value": "float64",
                "binom_bonferroni": "float64",
                "binom_fdr": "float64",
                "hyper_p_value": "float64",
                "hyper_bonferroni": "float64",
                "hyper_fdr": "float64",
            },
        )
        great_webserver.rename(columns={"ID": "id"}, inplace=True)
        del great_webserver["ontologie"]
        del great_webserver["term_name"]

        pp["before_pp_greatpy_size"].append(enrichment_tot.shape[0])
        enrichment_tot = enrichment_tot[enrichment_tot.index.isin(list(great_webserver["id"]))]
        pp["final_size"].append(enrichment_tot.shape[0])

        pp["before_pp_great_size"].append(great_webserver.shape[0])
        pp["%_of_GO_from_great_lost"].append(
            round(((great_webserver.shape[0] - enrichment_tot.shape[0]) / great_webserver.shape[0]) * 100, 2)
        )
        great_webserver = great_webserver[great_webserver["id"].isin(list(enrichment_tot.index))]

        great_webserver = great_webserver.sort_values("id")

        if disp_scatterplot or stats:
            binom_greatpy = []
            hyper_greatpy = []
            binom_great = []
            hyper_great = []
            for i in range(enrichment_tot.shape[0]):
                go_id = list(enrichment_tot.index)[i]
                curr_enrichment = enrichment_tot.iloc[i]
                curr_great_webserver = great_webserver.loc[great_webserver["id"] == go_id]
                binom_greatpy.append(float(curr_enrichment["binom_p_value"]))
                hyper_greatpy.append(float(curr_enrichment["hypergeom_p_value"]))
                binom_great.append(float(curr_great_webserver["binom_p_value"]))
                hyper_great.append(float(curr_great_webserver["hyper_p_value"]))
            binom = pd.DataFrame({"binom_greatpy": binom_greatpy, "binom_great": binom_great})
            hyper = pd.DataFrame({"hyper_greatpy": hyper_greatpy, "hyper_great": hyper_great})

            if disp_scatterplot:
                fig = plt.figure(figsize=(10, 5), dpi=80)
                fig.subplots_adjust(hspace=0.4, wspace=0.4)
                ax = fig.add_subplot(2, 2, 1)
                great.pl.scatterplot(binom, colname_x="binom_greatpy", colname_y="binom_great", title=None, ax=ax)
                ax = fig.add_subplot(2, 2, 2)
                great.pl.scatterplot(hyper, colname_x="hyper_greatpy", colname_y="hyper_great", title=None, ax=ax)
                fig.suptitle(f"results for {name}", fontsize=18)
                plt.show()

            if stats:
                stat_df["name"].append(name)
                stat_df["pearson_binom"].append(pearsonr(binom_great, binom_greatpy)[0])
                stat_df["pearson_hypergeom"].append(pearsonr(hyper_great, hyper_greatpy)[0])

        if good_gene_associations:
            gene_asso_great = pd.read_csv(
                great_asso,
                sep="\t",
                comment="#",
                names=["ontologies", "gene"],
                index_col=False,
                dtype={"ontologies": "object", "gene": "object"},
                usecols=["gene"],
            )
            gene_asso_greatpy = great.tl.get_association(
                test=pd.read_csv(
                    test,
                    sep="\t",
                    comment="#",
                    usecols=[0, 1, 2],
                    names=["chr", "chr_start", "chr_end"],
                    dtype={"chr": "object", "chr_start": "int64", "chr_end": "int64"},
                ),
                regdom=pd.read_csv(
                    regdom,
                    sep="\t",
                    comment="#",
                    names=["chr", "chr_start", "chr_end", "name", "tss", "strand"],
                    dtype={
                        "chr": "object",
                        "chr_start": "int64",
                        "chr_end": "int64",
                        "name": "object",
                        "tss": "int64",
                        "strand": "object",
                    },
                ),
            )

            in_in = gene_asso_great[gene_asso_great["gene"].isin(gene_asso_greatpy)].shape[0]
            in_out = [i for i in list(gene_asso_great["gene"]) if i not in gene_asso_greatpy]
            out_in = [i for i in gene_asso_greatpy if i not in list(gene_asso_great["gene"])]

            asso["name"].append(name)
            asso["number_good_gene_asso"].append(str(in_in))
            asso["number_genes_asso_lost"].append(str(len(in_out)))
            asso["number_gene_asso_excess"].append(str(len(out_in)))

    if stats and good_gene_associations:
        return pd.DataFrame(pp), pd.DataFrame(asso), pd.DataFrame(stat_df)

    elif stat_df:
        return pd.DataFrame(pp), pd.DataFrame(stat_df)

    elif good_gene_associations:
        pd.DataFrame(pp), pd.DataFrame(asso)
    else:
        return pd.DataFrame(pp)


def scale_data_5_75(data):
    mind = np.min(data)
    maxd = np.max(data)

    if maxd == mind:
        maxd = maxd + 1
        mind = mind - 1

    drange = maxd - mind
    return (((data - mind) / drange * 0.70) + 0.05) * 100


def plot_enrich(data: pd.DataFrame, n_terms: int = 20, color: str = "cool", save: bool = False):
    """
    Creation of a dotplot of the enrichment
    GO term in the inputs datas

    Parameters
    ----------
    data : pd.DataFrame
        Results for greatpy\n
    n_terms : int
        the number of term to be shown. \n
        Default is `20`\n
    color : str
        The color of the cmap in the plot\n
        Default is `"cool"`\n
    save : bool
        Is the plot should be save\n
        Default is `False`\n

    Returns
    -------
    None
        Dotplot of the enrichment n first GO term in the input data

    """
    # Test data input
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Please input a Pandas Dataframe output by gprofiler.")

    if not np.all([term in data.columns for term in ["p_value", "name", "intersection_size"]]):
        raise TypeError(f"The data frame {data} does not contain enrichment results from gprofiler.")

    data_to_plot = data.iloc[:n_terms, :].copy()
    data_to_plot["go.id"] = data_to_plot.index

    min_pval = data_to_plot["p_value"].min()
    max_pval = data_to_plot["p_value"].max()

    # Scale intersection_size to be between 5 and 75 for plotting
    # Note: this is done as calibration was done for values between 5 and 75
    data_to_plot["scaled.overlap"] = scale_data_5_75(data_to_plot["intersection_size"])

    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap=color, norm=norm)
    sm.set_array([])

    rcParams.update({"font.size": 14, "font.weight": "normal"})

    sns.set(style="whitegrid")

    path = plt.scatter(
        x="recall",
        y="name",
        c="p_value",
        cmap=color,
        norm=colors.LogNorm(min_pval, max_pval),
        data=data_to_plot,
        linewidth=1,
        edgecolor="grey",
        s=[(i + 10) ** 1.5 for i in data_to_plot["scaled.overlap"]],
    )
    ax = plt.gca()
    ax.invert_yaxis()

    ax.set_ylabel("")
    ax.set_xlabel("Gene ratio", fontsize=14, fontweight="normal")
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Get tick marks for this plot
    # Note: 6 ticks maximum
    min_tick = np.floor(np.log10(min_pval)).astype(int)
    max_tick = np.ceil(np.log10(max_pval)).astype(int)
    tick_step = np.ceil((max_tick - min_tick) / 6).astype(int)

    # Ensure no 0 values
    if tick_step == 0:
        tick_step = 1
        min_tick = max_tick - 1

    ticks_vals = [10**i for i in range(max_tick, min_tick - 1, -tick_step)]
    ticks_labs = ["$10^{" + str(i) + "}$" for i in range(max_tick, min_tick - 1, -tick_step)]

    # Colorbar
    fig = plt.gcf()
    cbaxes = fig.add_axes([0.8, 0.15, 0.03, 0.4])
    cbar = ax.figure.colorbar(sm, ticks=ticks_vals, shrink=0.5, anchor=(0, 0.1), cax=cbaxes)
    cbar.ax.set_yticklabels(ticks_labs)
    cbar.set_label("Adjusted p-value", fontsize=14, fontweight="normal")

    # Size legend
    min_olap = data_to_plot["intersection_size"].min()
    max_olap = data_to_plot["intersection_size"].max()
    olap_range = max_olap - min_olap

    # Note: approximate scaled 5, 25, 50, 75 values are calculated
    #      and then rounded to nearest number divisible by 5
    size_leg_vals = [
        np.ceil(i / 5) * 5
        for i in [min_olap, min_olap + (20 / 70) * olap_range, min_olap + (45 / 70) * olap_range, max_olap]
    ]
    size_leg_scaled_vals = scale_data_5_75(size_leg_vals)

    l1 = plt.scatter([], [], s=(size_leg_scaled_vals[0] + 10) ** 1.5, edgecolors="none", color="black")
    l2 = plt.scatter([], [], s=(size_leg_scaled_vals[1] + 10) ** 1.5, edgecolors="none", color="black")
    l3 = plt.scatter([], [], s=(size_leg_scaled_vals[2] + 10) ** 1.5, edgecolors="none", color="black")
    l4 = plt.scatter([], [], s=(size_leg_scaled_vals[3] + 10) ** 1.5, edgecolors="none", color="black")

    labels = [str(int(i)) for i in size_leg_vals]

    leg = plt.legend(
        [l1, l2, l3, l4],
        labels,
        ncol=1,
        frameon=False,
        fontsize=12,
        handlelength=1,
        loc="center left",
        borderpad=1,
        labelspacing=1.4,
        handletextpad=2,
        title="Gene overlap",
        scatterpoints=1,
        bbox_to_anchor=(-2, 1.5),
        facecolor="black",
    )

    if save:
        plt.savefig("dotplot_save.png", dpi=500)

    plt.show()


def make_bubble_heatmap(
    p_val_df: pd.DataFrame,
    odd_ratio_df: pd.DataFrame,
    na_color: str = "gray",
    tickscolorbar: list = [-2, -1, 0, 1, 2],
    vmin: int or float = -2.5,
    vmax: int or float = 2.5,
    heatmap_grid: list = [2, 4, 0, 2, 2, 1],
    circle_legend_grid: list = [2, 4, 0, 2, 2, 1],
    colorbar_grid: list = [2, 5, 0, 3, 2, 1],
    palette_id: str = "RdBu_r",
    cbar_label: str = "cbar_label",
    ncols: int or None = 8,
    marker: str or None = None,
    **kwargs,
) -> None:
    """
    Generate a dotplot with multiple categories

    Parameters
    ----------
    p_val_df : pd.DataFrame
        DataFrame of the p-values computed by enrichment function\n
    odd_ratio_df : pd.DataFrame
        DataFrame of the odds ratios computed by enrichment function\n
    na_color : str
        Color for NA values\n
        Default is `"gray"`\n
    tickscolorbar : list
        List of ticks should be add on the colorbar\n
        Default is `[-2, -1, 0, 1, 2]`\n
    vmin : int or float
        Minimum value for the colorbar\n
        Default is `-2.5`\n
    vmax : int or float
        Maximum value for the colorbar\n
        Default is `2.5`\n
    heatmap_grid : list
        Grid for the heatmap\n
        Default is `[2, 4, 0, 2, 2, 1]`\n
    circle_legend_grid : list
        Grid for the circle legend\n
        Default is `[2, 4, 0, 2, 2, 1]`\n
    colorbar_grid : list
        Grid for the colorbar\n
        Default is `[2, 5, 0, 3, 2, 1]`\n
    palette_id : str
        Palette id for the colorbar\n
        Default is `"RdBu_r"`\n
    cbar_label : str
        Label for the colorbar\n
        Default is `"cbar_label"`\n
    ncols : int or None
        Number of columns for the colorbar\n
        Default is `8`\n
    marker : str or None
        Marker for the dots\n
        Default is `None`\n
    kwargs
        Additional keyword arguments for the plot

    Returns
    -------
    None
        Bubble heatmap plot for gene ontology
    """
    from matplotlib import gridspec

    sns.set_style("white")

    rc("font", weight="normal")
    # rc("alpha",alpha=0.8)

    quantAmplifier = kwargs.get("quantAmplifier", 1)  # factor for size of bubbles

    # function you received on the color gradient
    pallete = None
    if palette_id == "RdBu_r":
        pallete = plt.cm.RdBu_r
    elif palette_id == "Blues":
        pallete = plt.cm.Blues
    elif palette_id == "Reds":
        pallete = plt.cm.Reds
    elif palette_id == "Greens":
        pallete = plt.cm.Greens
    elif palette_id == "Purples":
        pallete = plt.cm.Purples
    elif palette_id == "YlGnBu":
        pallete = plt.cm.YlGnBu
    elif palette_id == "PuOr":
        pallete = plt.cm.PuOr
    else:
        pass
    scalarmap, colorList = get_specific_color_gradient(pallete, np.array(p_val_df), vmin=vmin, vmax=vmax)

    fig = kwargs.get("fig", None)
    if fig is None:
        plt.clf()
        fig = plt.figure(figsize=(kwargs.get("w", 11), kwargs.get("h", 5)))

    nrows, ncols, rowi, coli, rowspan, colspan = heatmap_grid
    gs = gridspec.GridSpec(nrows, ncols)
    ax = plt.subplot(gs[rowi: rowi + rowspan, coli: coli + colspan])

    ylabelList = []
    i = 0

    # to keep the order of dataframes in the same order as sizes, revert the dataframes
    p_val_df = p_val_df.reindex(index=p_val_df.index[::-1])
    odd_ratio_df = odd_ratio_df.reindex(index=odd_ratio_df.index[::-1])

    min_circle_size = kwargs.get("min_circle_size")
    max_circle_size = kwargs.get("max_circle_size")

    for ri, r in p_val_df.iterrows():
        patList = list(r.values)

        sizeList = []
        if min_circle_size is not None and max_circle_size is not None:
            for si in list(odd_ratio_df.loc[ri]):
                if si < min_circle_size:
                    si = min_circle_size
                scaled = (si - min_circle_size) / (max_circle_size - min_circle_size)
                sizeList.append(scaled * quantAmplifier)
        else:
            sizeList = [(si ** kwargs.get("power", 1.0)) * (quantAmplifier) for si in list(odd_ratio_df.loc[ri])]

        colorList = scalarmap.to_rgba(patList)
        colorList = [ci if not np.isnan(vi) else na_color for ci, vi in zip(colorList, patList)]
        x = list(range(len(patList)))
        y = [i] * len(patList)

        # define hataches, linewidths and alphas
        hatches = [None if np.isnan(pi) else None for pi in patList]
        linewidths = (
            [0.1 if np.isnan(pi) else kwargs.get("sig_line_width", 2.0) for pi in patList]
            if kwargs.get("line_widths", None) is None
            else list(kwargs.get("line_widths").loc[ri])
        )
        alphas = [0.5 if np.isnan(alpha) else 1.0 for alpha in patList]

        for xi, yi, si, ci, hatch_i, lw, alpha in zip(x, y, sizeList, colorList, hatches, linewidths, alphas):
            ax.scatter(
                xi,
                yi,
                s=abs(si) * quantAmplifier,
                marker=marker if marker is not None else ("v" if si < 0 else "^"),
                hatch=hatch_i,
                alpha=alpha,
                edgecolor=kwargs.get("edgecolor", "black"),
                color=ci,
                linewidth=lw,
            )
        ylabelList.append(ri)
        i += 1
    if kwargs.get("grid", True):
        plt.grid(True, linewidth=kwargs.get("grid_linewidth", 1.0))
    plt.title(kwargs.get("heatmap_title", "title"), fontname="Arial", fontweight="normal")
    plt.xlabel(kwargs.get("xlab", "xlab"), fontsize=14, fontname="Arial", fontweight="normal")
    plt.ylabel(kwargs.get("ylab", "ylab"), fontsize=14, fontname="Arial", fontweight="normal")

    plt.yticks(list(range(len(ylabelList))))
    ax.set_yticklabels(
        ylabelList,
        fontsize=kwargs.get("yticks_fontsize", 11),
        color="black",
        ha=kwargs.get("ha_ylabs", "right"),
        fontname="Arial",
        fontweight="normal",
        alpha=1,
    )

    remove_top_n_right_ticks(ax)

    plt.xticks(
        list(range(len(p_val_df.columns))),
        p_val_df.columns,
        fontsize=kwargs.get("xticks_fontsize", 11),
        rotation=kwargs.get("rotation_xlabs", 90),
        color="black",
        ha=kwargs.get("ha_xlabs", "center"),
        fontname="Arial",
        fontweight="normal",
    )

    lh, lt = get_legendHandle_for_second_sanity_check_plot(
        quantAmplifier=quantAmplifier * quantAmplifier,
        marker="^" if marker is None else marker,
        fmt=kwargs.get("fmt_legend", "%.1f"),
        lw=0.5,
        min_circle_size=min_circle_size,
        max_circle_size=max_circle_size,
        values=[tick ** kwargs.get("power", 1.0) for tick in kwargs.get("circle_legend_ticks")],
        labels=kwargs.get("circle_legend_ticks"),
    )
    if kwargs.get("show_circle_legend", True):
        l = plt.legend(
            lh,
            lt,
            bbox_to_anchor=kwargs.get("circle_legend_bbox", (1.8, 1)),
            scatterpoints=1,
            title=kwargs.get("circles_legend_title", "circles_legend_title"),
            ncol=1,
            frameon=False,
        )
        # Add the legend manually to the current Axes.
        ax = plt.gca().add_artist(l)

    # this is to add the circle (significant or not
    if kwargs.get("show_sig_legend", False):
        lh, lt = get_legendHandle_for_second_sanity_check_plot(
            quantAmplifier=quantAmplifier * quantAmplifier,
            marker="^" if marker is None else marker,
            fmt=kwargs.get("sig_fmt_legend", "%s"),
            lw=[kwargs.get("sig_line_width", 2.0), 0.0],
            labels=kwargs.get("sig_legend_ticks", ["Yes", "No"]),
            values=kwargs.get("sig_legend_ticks", ["Yes", "No"]),
            edgecolor=kwargs.get("edgecolor", "black"),
            min_size_default=max(kwargs.get("circle_legend_ticks")) ** kwargs.get("power", 1.0),
        )

        plt.legend(
            lh,
            lt,
            bbox_to_anchor=kwargs.get("sig_legend_bbox", (2.4, 1)),
            title=kwargs.get("sig_legend_title", "significant"),
            ncol=1,
            scatterpoints=1,
            frameon=False,
        )

    if kwargs.get("show_colorbar", True):
        cbar = fig.colorbar(
            scalarmap,
            orientation="horizontal",
            format=kwargs.get("cbar_fmt_ticks", "%.1f"),
            ticks=tickscolorbar,
            pad=0.1,
        )
        cbar.ax.tick_params(labelsize=kwargs.get("colorbar_ticks_labelsize", 12))
        cbar.set_label(cbar_label, fontsize=12, fontweight="normal")
    despine_all()


def get_specific_color_gradient(colormap, inputList, **kwargs) -> tuple:
    vmin = kwargs.get("vmin", "blaq")
    vmax = kwargs.get("vmax", "blaq")
    cm = plt.get_cmap(colormap)
    if vmin == "blaq" or vmax == "blaq":
        if type(inputList) == list:
            cNorm = matplotlib.colors.Normalize(vmin=min(inputList), vmax=max(inputList))
        else:
            cNorm = matplotlib.colors.Normalize(vmin=inputList.min(), vmax=inputList.max())
    else:
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    scalarMap.set_array(inputList)
    colorList = scalarMap.to_rgba(inputList)
    return scalarMap, colorList


def despine_all() -> None:
    sns.despine(offset=10, trim=True, top=True, right=True, left=True, bottom=True)


def remove_top_n_right_ticks(ax) -> None:
    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")


def get_legendHandle_for_second_sanity_check_plot(
    quantAmplifier=None,
    labels=None,
    values=None,
    marker="^",
    fmt="%.1f",
    lw=2.0,
    min_circle_size=None,
    max_circle_size=None,
    min_size_default=1.0,
    edgecolor="black",
    color="grey",
):
    labels = [0.1, 0.5, 1.0, 1.5] if labels is None else labels

    legendHandleList = list()
    labelSize = list()
    for i, quantMem in enumerate(values):
        values[i]

        if min_circle_size is not None and max_circle_size is not None:
            scaled = (quantMem - min_circle_size) / (max_circle_size - min_circle_size)
            quantMem = scaled

        next_size = (quantMem * quantAmplifier) if not isinstance(quantMem, str) else min_size_default * quantAmplifier

        legendHandleList.append(
            plt.scatter(
                [],
                [],
                s=next_size,
                color=color,
                edgecolor=edgecolor,
                linewidths=lw if not isinstance(lw, list) else lw[i],
                alpha=0.9,
                marker=marker,
            )
        )
        labelSize.append(fmt % labels[i])
    return legendHandleList, labelSize


def dotplot_multi_sample(
    test_data: dict,
    n_row: int = 5,
    list_id: list = [],
    fig: matplotlib.figure.Figure or None = None,
    show_term_name: bool = False,
    term_name_nchars: int = 30,
    dot_size_amplifier: int or float = 7,
    palette_id: str or matplotlib.cm = "Reds",
    ylab: str = "GO",
    xlab: str = "",
    label_colorbar: str = "-log(p_hypergeometric)",
    marker: str = "o",
    plot_title: str = "Dotplot of enrichment GO terms",
    line_width: int = 0.1,
    circle_legend: str = "log2(odd ratio)",
    save: bool = False,
    **kwargs,
):
    """
    Dotplot of enrichment GO terms for a given list of example genomic regions.

    Parameters
    ----------
    test_data : dict
        dict of multiple tests output\n
    n_row : int
        Number of rows pick in each dataframe.\n
        Default is `5`\n
    list_id : list
        List of IDs to be plotted. If `list_id == []`, any filter will be applied.\n
        Default is `[]`\n
    fig : matplotlib.figure.Figure or None
        Figure to plot the dotplot.\n
        Default is `None`\n
    show_term_name : bool
        Whether to show the GO term name.\n
        Default is `False`\n
    term_name_nchars : int
        Number of characters to show for the GO term name.\n
        Default is `30`\n
    dot_size_amplifier : int or float
        Amplifier for the dot size.\n
        Default is `7`\n
    palette_id : str or matplotlib.cm
        Color palette for the dotplot.\n
        Default is `Reds`\n
    ylab : str
        Y-axis label.\n
        Default is `"GO"`\n
    xlab : str
        X-axis label.\n
        Default is `""`\n
    label_colorbar : str
        Label for the colorbar.\n
        Default is `"-log(p_hypergeometric)"`\n
    marker : str
        Marker for the dotplot.\n
        Default is `"o"`\n
    plot_title : str
        Plot title.\n
        Default is `"Dotplot of enrichment GO terms"`\n
    line_width : int
        Dot line width.\n
        Default is `0.1`\n
    circle_legend : str
        Legend for the circle.\n
        Default is `"log2(odd ratio)"`\n
    kwargs
        Other parameters to be passed to the make make_bubble_heatmap function.

    Returns
    -------
    None
        Dotplot of enrichment GO terms for the given df\n
    p_val : pandas.DataFrame
        Dataframe of plotted p-values\n
    odds_ratio : pandas.DataFrame
        Dataframe of plotted odds ratios\n
    df : pandas.DataFrame
        Dataframe of all results concatenated\n

    Examples
    --------
    >>> test = ["SRF:Ishikawa,A-673-clone-Asp114,K-562,MCF-7,Hep-G2","MAX:K-562,WA01,HeLa-S3", "BACH1:A-549,GM12878"]
    >>> tmp_df = great.tl.enrichment_multiple(
        tests = test, regdom_file="../data/human/hg38/regulatory_domain.bed",
        chr_size_file="../data/human/hg38/chr_size.bed",
        annotation_file="../data/human/ontologies.csv", binom=True, hypergeom=True,)
    >>> p_val,odd_ratio = great.pl.dotplot_multi_sample(tmp_df)

    >>> p_val
    ...    | id         |       0 |       1 |        2 |
    ...    |:-----------|--------:|--------:|---------:|
    ...    | GO:0051292 | 6.28405 | 1       |  1       |
    ...    | GO:0030261 | 5.51863 | 1       |  1       |
    ...    | GO:0001650 | 5.3019  | 1       |  1       |
    ...    | GO:0090096 | 5.29709 | 1       |  1       |
    ...    | GO:0099637 | 5.29709 | 1       |  1       |
    ...    | GO:0004896 | 1       | 7.98072 |  1       |
    ...    | GO:0038165 | 1       | 7.93887 |  1       |
    ...    | GO:0030883 | 1       | 7.03282 |  1       |
    ...    | GO:0048006 | 1       | 7.03282 |  1       |
    ...    | GO:0004924 | 1       | 7.03282 |  1       |
    ...    | GO:0008137 | 1       | 1       | 14.0412  |
    ...    | GO:0015990 | 1       | 1       | 13.2519  |
    ...    | GO:0006120 | 1       | 1       | 11.7472  |
    ...    | GO:0045277 | 1       | 1       | 10.2356  |
    ...    | GO:0030964 | 1       | 1       |  9.74806 |

    >>> odds_ratio
    ...    | id         |       0 |       1 |       2 |
    ...    |:-----------|--------:|--------:|--------:|
    ...    | GO:0051292 | 4.94165 | 0       | 0       |
    ...    | GO:0030261 | 4.39416 | 0       | 0       |
    ...    | GO:0001650 | 2.54406 | 0       | 0       |
    ...    | GO:0090096 | 7.64209 | 0       | 0       |
    ...    | GO:0099637 | 7.64209 | 0       | 0       |
    ...    | GO:0004896 | 0       | 3.58059 | 0       |
    ...    | GO:0038165 | 0       | 6.00685 | 0       |
    ...    | GO:0030883 | 0       | 5.42189 | 0       |
    ...    | GO:0048006 | 0       | 5.42189 | 0       |
    ...    | GO:0004924 | 0       | 5.42189 | 0       |
    ...    | GO:0008137 | 0       | 0       | 3.76802 |
    ...    | GO:0015990 | 0       | 0       | 6.61172 |
    ...    | GO:0006120 | 0       | 0       | 3.65752 |
    ...    | GO:0045277 | 0       | 0       | 5.44179 |
    ...    | GO:0030964 | 0       | 0       | 7.02675 |

    """
    min_p, min_od, max_p, max_od = 100, 100, 0, 0
    gene_tot = {}
    df = pd.DataFrame()

    for i, name in enumerate(test_data.keys()):
        gene_tot[i] = name
        res = test_data[name].sort_values(by="hypergeom_p_value")
        if len(list_id) > 0:
            res = res[res.index.isin(list_id)]
        res = res.iloc[:n_row]

        res = res.reset_index().rename(columns={"index": "id"})
        res["index"] = [i] * res.shape[0]
        frames = [df, res]
        df = pd.concat(frames)
    df["hypergeom_p_value"] = -np.log(df["hypergeom_p_value"])
    df["hypergeometric_fold_enrichment"] = np.log2(df["hypergeometric_fold_enrichment"])
    dup = df[df.duplicated("id") == True]
    duplicate = dup.drop(
        columns=["go_term", "binom_p_value", "binom_fold_enrichment", "intersection_size", "recall"]
    ).to_dict()
    df = df.drop_duplicates("id")
    p_val = df.pivot("id", "index", "hypergeom_p_value").fillna(1).reindex(df["id"])
    odd_ratio = df.pivot("id", "index", "hypergeometric_fold_enrichment").fillna(0).reindex(df["id"])

    if dup.shape[0] > 0:
        for i, name in enumerate(duplicate["id"].values()):
            if duplicate["index"][i] in p_val.columns:
                p_val[duplicate["index"][i]][name] = duplicate["hypergeom_p_value"][i]
                odd_ratio[duplicate["index"][i]][name] = duplicate["hypergeometric_fold_enrichment"][i]
            else:
                p_val[duplicate["index"][i]] = duplicate["hypergeom_p_value"][i]
                odd_ratio[duplicate["index"][i]] = duplicate["hypergeometric_fold_enrichment"][i]

    p_val.rename(columns=gene_tot, inplace=True)
    odd_ratio.rename(columns=gene_tot, inplace=True)

    if show_term_name:
        go_n = []
        for i in range(df.shape[0]):
            curr = df.iloc[i]
            go_n.append(curr["id"] + " " + curr["go_term"][:term_name_nchars])
        p_val["test"] = go_n
        odd_ratio["test"] = go_n

        p_val = p_val.set_index("test")
        odd_ratio = odd_ratio.set_index("test")

    for i in p_val.columns:
        if min_p > p_val[i].min():
            min_p = p_val[i].min()
        if max_p < p_val[i].max():
            max_p = p_val[i].max()

    for i in odd_ratio.columns:
        if min_od > odd_ratio[i].min():
            min_od = odd_ratio[i].min()
        if max_od < odd_ratio[i].max():
            max_od = odd_ratio[i].max()

    max_p = int(round(max_p))
    max_od = int(round(max_od))
    min_p = int(round(min_p))
    min_od = int(round(min_od))

    if min_od == 0:
        min_od = 1

    plt.rcParams.update({"font.size": 14, "font.weight": "normal"})
    plt.tight_layout(pad=0.1, h_pad=0.1, w_pad=0.1)

    make_bubble_heatmap(
        p_val,
        odd_ratio,
        quantAmplifier=dot_size_amplifier,
        circle_legend_ticks=[i for i in range(round(min_od), round(max_od) + 1, 2)],
        palette_id=palette_id,
        ylab=ylab,
        xlab=xlab,
        tickscolorbar=[int(i) for i in range(round(min_p), round(max_p) + 1, 2)],
        vmin=min_p,
        vmax=max_p,
        cbar_label=label_colorbar,
        marker=marker,
        heatmap_title=plot_title,
        sig_line_width=line_width,
        circles_legend_title=circle_legend,
        cbar_fmt_ticks="%.0f",
        fig=fig,
    )
    if save:
        plt.savefig("multidotplot.png", dpi=300)

    return p_val, odd_ratio, df
