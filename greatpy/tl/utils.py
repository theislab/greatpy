import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

import greatpy as gp
import greatpy as great


def get_nb_asso_per_region(test: str or pd.DataFrame, regdom: str or pd.DataFrame) -> dict:
    """
    From a file of genomic regions from CHIPseq
    and a file of genomic regulatory domains, determine number of peaks
    associated with each gene in the regulatory domain.

    Parameters
    ----------
    test : str or pd.DataFrame
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : str or pd.DataFrame
        path of the file with the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"Name"	"tss"	"strand"].

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
            "Name":["RNF223","C1orf159"],
            "tss":[1074306,1116089],
            "strand":['-','-']
        })
    >>> get_association(test,regdom)
    ...    {'RNF223':2}

    """
    res = {}
    test, regdom, _, _ = gp.tl.GREAT.loader(test, regdom, None, None)
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
    From a file of genomic regions from CHIPseq
    and a file of genomic regulatory domains, determine the distance from peaks
    to the transcription start site of the associated gene

    Parameters
    ----------
    test : str
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : str
        path of the file with the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"Name"	"tss"	"strand"].

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
            "Name":["RNF223","C1orf159"],
            "tss":[1074306,1116089],
            "strand":['-','-']
        })
    >>> get_association(test,regdom)
    ...    {'RNF223':[-22278]}

    """
    res = {}
    test, regdom, _, _ = gp.tl.GREAT.loader(test, regdom, None, None)
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


def get_all_comparison(
    preprocessing: bool = True, good_gene_associations: bool = True, disp_scatterplot: bool = True, stats: bool = True
):
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

    len(os.listdir("../data/tests/test_data/input/"))
    for path in os.listdir("../data/tests/test_data/input/"):
        sp = path.split(".")
        id = sp[0][:2]
        name = sp[0][3:]
        pp["name"].append(name)
        i = 0
        great_out = ""
        great_asso = ""

        for out_path in os.listdir("../data/tests/test_data/output/"):
            if out_path.split("_")[0] == id:
                if re.match(".*hg19.*", out_path) != None:
                    assembly = "hg19"
                else:
                    assembly = "hg38"

                if re.match(".*output.*", out_path) != None:
                    great_out = "../data/tests/test_data/output/" + out_path
                else:
                    great_asso = "../data/tests/test_data/output/" + out_path

        test = "../data/tests/test_data/input/" + path
        regdom = f"../data/human/{assembly}/regulatory_domain.bed"
        size = f"../data/human/{assembly}/chr_size.bed"

        if great_out == "" or great_asso == "":
            return False

        enrichment_tot = gp.tl.GREAT.enrichment(
            test_file=test,
            regdom_file=regdom,
            chr_size_file=size,
            annotation_file=f"../data/human/ontologies.csv",
            binom=True,
            hypergeom=True,
        )
        enrichment_tot = gp.tl.GREAT.set_bonferroni(enrichment_tot, 0.05)
        enrichment_tot = gp.tl.GREAT.set_fdr(enrichment_tot, 0.05)

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
                gp.pl.scatterplot(binom, colname_x="binom_greatpy", colname_y="binom_great", title=None, ax=ax)
                ax = fig.add_subplot(2, 2, 2)
                gp.pl.scatterplot(hyper, colname_x="hyper_greatpy", colname_y="hyper_great", title=None, ax=ax)
                fig.suptitle(f"results for {name}")
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
                    names=["Chr", "Chr_Start", "Chr_End"],
                    dtype={"Chr": "object", "Chr_Start": "int64", "Chr_End": "int64"},
                ),
                regdom=pd.read_csv(
                    regdom,
                    sep="\t",
                    comment="#",
                    names=["Chr", "Chr_Start", "Chr_End", "Name", "tss", "Strand"],
                    dtype={
                        "Chr": "object",
                        "Chr_Start": "int64",
                        "Chr_End": "int64",
                        "Name": "object",
                        "tss": "int64",
                        "Strand": "object",
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

    if stat_df:
        stat_df = pd.DataFrame(stat_df)

    return pd.DataFrame(pp), pd.DataFrame(asso), stat_df
