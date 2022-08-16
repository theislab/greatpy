import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2
import seaborn as sns
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as r
from rpy2.robjects.packages import importr
from scipy.stats import pearsonr

import greatpy as great

pandas2ri.activate()


def get_nb_asso_per_region(test: str or pd.DataFrame, regdom: str or pd.DataFrame) -> dict:
    """
    Determine number of peaks associated with each gene in the regulatory domain.

    Parameters
    ----------
    test : str or pd.DataFrame
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : str or pd.DataFrame
        path of the file with the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"name"	"tss"	"strand"].

    Returns
    -------
    res : dict
        dict with the number of associated genes per genomic region : key = associated gene, value = number of peaks associated with the gene

    Examples
    --------
    >>> test = pd.DataFrame(
        {
            "chr":["chr1"],
            "chr_start":[1052028],
            "chr_end": [1052049]}
        )
    >>> regdom = pd.DataFrame(
        {
            "chr":["chr1","chr1"],
            "chr_start":[1034992,1079306],
            "chr_end": [1115089,1132016],
            "name":["RNF223","C1orf159"],
            "tss":[1074306,1116089],
            "strand":['-','-']
        })
    >>> get_association(test,regdom)
    ...    {'RNF223':2}

    """
    res = {}
    test, regdom, _, _ = great.tl.GREAT.loader(test, regdom, None, None)
    for i in range(test.shape[0]):
        currTest = test.iloc[i]
        regdom_curr_test = regdom.loc[(regdom["chr"] == currTest["chr"])].sort_values("chr_start")
        regdom_curr_test = regdom_curr_test.loc[
            (
                (regdom_curr_test["chr_start"] <= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] >= currTest["chr_end"])
            )
            | (  # regdom overlap totally test
                (regdom_curr_test["chr_start"] >= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] <= currTest["chr_end"])
            )
            | (  # test overlap totally regdom
                (regdom_curr_test["chr_start"] <= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] <= currTest["chr_end"])
                & (regdom_curr_test["chr_end"] >= currTest["chr_start"])
            )
            | (  # regdom overlap not totally test on left side
                (regdom_curr_test["chr_start"] >= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] >= currTest["chr_end"])
                & (regdom_curr_test["chr_start"] <= currTest["chr_end"])
            )  # regdom overlap not totally test on right side
        ]
        res[i] = regdom_curr_test.shape[0]
    return res


def get_dist_to_tss(test: str or pd.DataFrame, regdom: str or pd.DataFrame) -> dict:
    """
    Determine the distance from peaks to the transcription start site of the associated gene

    Parameters
    ----------
    test : str or pd.DataFrame
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : str or pd.DataFrame
        path of the file with the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"name"	"tss"	"strand"].

    Returns
    -------
    res : dict
        dict with the distance from tss to the associated genes : key = number of the input, value = distance from peaks to tss of associated genes

    Examples
    --------
    >>> test = pd.DataFrame(
        {
            "chr":["chr1"],
            "chr_start":[1052028],
            "chr_end": [1052049]}
        )
    >>> regdom = pd.DataFrame(
        {
            "chr":["chr1","chr1"],
            "chr_start":[1034992,1079306],
            "chr_end": [1115089,1132016],
            "name":["RNF223","C1orf159"],
            "tss":[1074306,1116089],
            "strand":['-','-']
        })
    >>> get_association(test,regdom)
    ...    {'RNF223':[-22278]}

    """
    res = {}
    test, regdom, _, _ = great.tl.GREAT.loader(test, regdom, None, None)
    for i in range(test.shape[0]):
        currTest = test.iloc[i]
        mean_pos_test = (currTest["chr_end"] + currTest["chr_start"]) / 2
        regdom_curr_test = regdom.loc[(regdom["chr"] == currTest["chr"])].sort_values("chr_start")
        regdom_curr_test = regdom_curr_test.loc[
            (
                (regdom_curr_test["chr_start"] <= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] >= currTest["chr_end"])
            )
            | (  # regdom overlap totally test
                (regdom_curr_test["chr_start"] >= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] <= currTest["chr_end"])
            )
            | (  # test overlap totally regdom
                (regdom_curr_test["chr_start"] <= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] <= currTest["chr_end"])
                & (regdom_curr_test["chr_end"] >= currTest["chr_start"])
            )
            | (  # regdom overlap not totally test on left side
                (regdom_curr_test["chr_start"] >= currTest["chr_start"])
                & (regdom_curr_test["chr_end"] >= currTest["chr_end"])
                & (regdom_curr_test["chr_start"] <= currTest["chr_end"])
            )  # regdom overlap not totally test on right side
        ]
        res[i] = []
        for j in range(regdom_curr_test.shape[0]):
            res[i].append(regdom_curr_test.iloc[j]["tss"] - mean_pos_test)
    return res


def online_vs_local_vs_greatpy_comparison(
    input_folder:str = "../data/tests/test_data/input/", 
    information_folder:str = "../data/human/",
    annotation_file: str = "../../data/human/ontologies.csv"
    ) :
    """
    Make a comparison between the online and the local version of rGREAT and greatpy.
    The function return a clustermap of the results between online vs local and greatpy.

    Parameters
    ----------
    input_folder : str
        path of the folder with the input files for the tests.\n
        Default is `"../data/tests/test_data/input/"`. \n
    information_folder : str
        path of the folder with the information files for the tests. \n
        Default is `"../data/human/"`\n
        The input folder should contains the files : \n
        - information_folder/assembly_eg_hg38/regulatory_domain.bed \n
        - information_folder/assembly_eg_hg38/chr_size.bed\n
    annotation_file : str
        path of the file with the ontologies annotations. The file should be in the csv format with `;` as separator and contain the following columns :\n
        - "ensembl": ensemble id (optionnal)\n
        - "id": id of the ontology\n
        - "name": name of the ontology\n
        - "ontology.group": group of the ontology\n
        - "gene.name": name of the gene\n
        - "symbol": symbol of the gene equivalent to the gene.name\n
        Default is `"../../data/human/ontologies.csv"`\n

    Returns
    -------
    None

    Note
    ----
    To use this function, you should have installed in your environment: 
        - rpy2
        - R base with the version 3.6.1 
        - The following R packages : `rGREAT`, `GenomicRanges`
    
    """
    importr("rGREAT")
    ranges = importr("GenomicRanges")

    stat_df = {"name": [], "pearson_binom": [], "pearson_hypergeom": []}
    pp = {
        "name": [],
        "before_pp_greatpy_size": [],
        "before_pp_local_size": [],
        "final_size": [],
        "%_of_diffrent_GO_term": [],
    }

    ann = pd.read_csv(
        annotation_file,
        sep=";",
        comment="#",
        header=0,
    )
    ann["name"] = ann["name"].str.lower()
    pd.DataFrame()

    res = {"name1": [], "file": [], "name3": [], "value": []}

    for name in [
        "01_random.bed",
        "04_ultra_hg38.bed",
        "06_height_snps_hg38.bed",
        "07_height_snps_hg19.bed",
        "10_MAX.bed",
    ]:
        # find the assembly
        if re.match(".*hg19.*", name) != None:
            assembly = "hg19"
        else:
            assembly = "hg38"

        # online test
        res_online = rpy2.robjects.r["submitGreatJob"](
            f"{input_folder}{name}", species=f"{assembly}", help=False
        )
        res_online = rpy2.robjects.r["getEnrichmentTables"](res_online)

        # local test
        # proprocessing : make a Grange frame
        df = r["read.csv"](f"{input_folder}{name}", sep="\t")
        seqname = rpy2.robjects.StrVector(
            ["seqnames", "seqname", "chromosome", "X.Chr", "chr", "chromosome_name", "seqid"]
        )
        rpy2.robjects.StrVector(["end", "stop"])
        df = ranges.makeGRangesFromDataFrame(df, seqnames_field=seqname)

        # great calculation
        local = rpy2.robjects.r["great"](df, "GREAT:C5", f"txdb:{assembly}", verbose=False)
        local = rpy2.robjects.r["getEnrichmentTables"](local)

        # greatpy calculation
        greatpy = great.tl.GREAT.enrichment(
            test_file=f"{input_folder}{name}",
            regdom_file=f"{information_folder}{assembly}/regulatory_domain.bed",
            chr_size_file=f"{information_folder}{assembly}/chr_size.bed",
            annotation_file=f"{annotation_file}",
            binom=True,
            hypergeom=True,
        )

        # create each dataframe
        # online
        name_online = [
            cdc.lower()
            for cdc in list(res_online.rx2("GO Molecular Function").rx2("name"))
            + list(res_online.rx2("GO Biological Process").rx2("name"))
            + list(res_online.rx2("GO Cellular Component").rx2("name"))
        ]
        online = pd.DataFrame(
            {
                "id": list(res_online.rx2("GO Molecular Function").rx2("ID"))
                + list(res_online.rx2("GO Biological Process").rx2("ID"))
                + list(res_online.rx2("GO Cellular Component").rx2("ID")),
                "name": name_online,
                "binom_p_val": list(res_online.rx2("GO Molecular Function").rx2("Binom_Raw_PValue"))
                + list(res_online.rx2("GO Biological Process").rx2("Binom_Raw_PValue"))
                + list(res_online.rx2("GO Cellular Component").rx2("Binom_Raw_PValue")),
                "hyper_p_val": list(res_online.rx2("GO Molecular Function").rx2("Hyper_Raw_PValue"))
                + list(res_online.rx2("GO Biological Process").rx2("Hyper_Raw_PValue"))
                + list(res_online.rx2("GO Cellular Component").rx2("Hyper_Raw_PValue")),
            }
        )

        # local
        name_local = list(local.rx2("id"))
        name_local = [" ".join(cdc.lower().split("_")[1:]) for cdc in list(local.rx2("id"))]
        local = pd.DataFrame(
            {
                "name": name_local,
                "binom_p_val": list(local.rx2("p_value")),
                "hyper_p_val": list(local.rx2("p_value_hyper")),
            }
        )

        # greatpy
        greatpy["go_term"] = greatpy["go_term"].str.lower()

        # correlation between online and greatpy
        o_great = online[online["id"].isin(list(greatpy.index))]
        greatpy[greatpy.index.isin(list(online["id"]))]

        tot = pd.DataFrame()
        res_bin = []
        res_hyp = []
        cols_to_supr = []
        for i in range(len(o_great)):
            curr_o = o_great.iloc[i]
            if greatpy.loc[greatpy.index == curr_o["id"]].shape[0] > 0:
                res_bin.append(float(greatpy.loc[greatpy.index == curr_o["id"]]["binom_p_value"]))
                res_hyp.append(float(greatpy.loc[greatpy.index == curr_o["id"]]["hypergeom_p_value"]))
            else:
                cols_to_supr.append(i)

        tot["id"] = o_great["id"]
        tot["binom_p_value_online"] = -np.log(list(o_great["binom_p_val"]))
        tot["hyper_p_value_online"] = -np.log(list(o_great["hyper_p_val"]))
        tot["binom_p_value_greatpy"] = -np.log(res_bin)
        tot["hypergeom_p_value_greatpy"] = -np.log(res_hyp)

        # add GO id to local
        go = []
        for i in range(local.shape[0]):
            curr = local.iloc[i]

            if ann.loc[ann["name"].isin([curr["name"]])].shape[0] > 0:
                go.append(ann.loc[ann["name"].isin([curr["name"]])].iloc[0]["id"])
            else:
                go.append("")
        local["id"] = go
        local = local.loc[local["id"] != ""]

        # Correlation between online and local
        o_loc = online[online["id"].isin(local["id"])]
        loc_o = local[local["id"].isin(online["id"])]
        tot_2 = pd.DataFrame()

        res_bin = []
        res_hyp = []
        cols_to_supr = []
        for i in range(len(o_loc)):
            curr_o = o_loc.iloc[i]
            if loc_o.loc[loc_o["id"] == curr_o["id"]].shape[0] > 0:
                res_bin.append(float(loc_o.loc[loc_o["id"] == curr_o["id"]].iloc[0]["binom_p_val"]))
                res_hyp.append(float(loc_o.loc[loc_o["id"] == curr_o["id"]].iloc[0]["hyper_p_val"]))
            else:
                cols_to_supr.append(i)
        tot_2["id"] = o_loc["id"]
        tot_2["binom_p_value_online"] = -np.log(list(o_loc["binom_p_val"]))
        tot_2["hyper_p_value_online"] = -np.log(list(o_loc["hyper_p_val"]))
        tot_2["binom_p_value_local"] = -np.log(res_bin)
        tot_2["hypergeom_p_value_local"] = -np.log(res_hyp)

        # Correlation between greatpy and local
        g_loc = greatpy[greatpy.index.isin(local["id"])]
        loc_g = local[local["id"].isin(greatpy.index)]
        tot_3 = pd.DataFrame()

        res_bin = []
        res_hyp = []
        cols_to_supr = []
        g_loc = g_loc.reset_index().rename(columns={"index": "id"})
        for i in range(g_loc.shape[0]):
            curr_o = g_loc.iloc[i]
            if loc_g.loc[loc_g["id"] == curr_o["id"]].shape[0] > 0:
                res_bin.append(float(loc_g.loc[loc_g["id"] == curr_o["id"]].iloc[0]["binom_p_val"]))
                res_hyp.append(float(loc_g.loc[loc_g["id"] == curr_o["id"]].iloc[0]["hyper_p_val"]))
            else:
                cols_to_supr.append(i)
        tot_3["id"] = g_loc["id"]
        tot_3["binom_p_value_greatpy"] = -np.log(list(g_loc["binom_p_value"]))
        tot_3["hyper_p_value_greatpy"] = -np.log(list(g_loc["hypergeom_p_value"]))
        tot_3["binom_p_value_local"] = -np.log(res_bin)
        tot_3["hypergeom_p_value_local"] = -np.log(res_hyp)

        tot = tot.replace(np.inf, np.nan)
        tot_2 = tot_2.replace(np.inf, np.nan)
        tot_3 = tot_3.replace(np.inf, np.nan)

        tot = tot.dropna()
        tot_2 = tot_2.dropna()
        tot_3 = tot_3.dropna()

        # create the table
        name_plot = " ".join(name.split(".")[0].split("_")[1:])
        res["name1"].append("online")
        res["file"].append(name_plot)
        res["name3"].append("greatpy")
        res["value"].append(pearsonr(tot["binom_p_value_online"], tot["binom_p_value_greatpy"])[0])

        res["name1"].append("online")
        res["file"].append(name_plot)
        res["name3"].append("local")
        res["value"].append(pearsonr(tot_2["binom_p_value_online"], tot_2["binom_p_value_local"])[0])

    res = pd.DataFrame(res)
    res = res.pivot("file", "name3", "value")
    g = sns.heatmap(data=res, cmap="Reds", annot=True)
    g.set_title("correlation with GREAT server")
    g.set_ylabel("online results for each file")
    g.set_xlabel("algorithm results for each file")

    plt.show(g)
