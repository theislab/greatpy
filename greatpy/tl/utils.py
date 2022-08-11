import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

import greatpy as great


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
    test : str
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : str
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


def get_all_comparison(good_gene_associations: bool = True, disp_scatterplot: bool = True, stats: bool = True):
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

        enrichment_tot = great.tl.GREAT.enrichment(
            test_file=test,
            regdom_file=regdom,
            chr_size_file=size,
            annotation_file="../data/human/ontologies.csv",
            binom=True,
            hypergeom=True,
        )
        enrichment_tot = great.tl.GREAT.set_bonferroni(enrichment_tot, 0.05)
        enrichment_tot = great.tl.GREAT.set_fdr(enrichment_tot, 0.05)

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
                    names=["Chr", "chr_start", "chr_end"],
                    dtype={"Chr": "object", "chr_start": "int64", "chr_end": "int64"},
                ),
                regdom=pd.read_csv(
                    regdom,
                    sep="\t",
                    comment="#",
                    names=["Chr", "chr_start", "chr_end", "name", "tss", "Strand"],
                    dtype={
                        "Chr": "object",
                        "chr_start": "int64",
                        "chr_end": "int64",
                        "name": "object",
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

    if stat_df and good_gene_associations:
        return pd.DataFrame(pp), pd.DataFrame(asso), pd.DataFrame(stat_df)

    elif stat_df:
        return pd.DataFrame(pp), pd.DataFrame(stat_df)

    elif good_gene_associations:
        pd.DataFrame(pp), pd.DataFrame(asso)
    else:
        return pd.DataFrame(pp)


def online_vs_local_vs_greatpy_comparison():
    stat_df = {"name": [], "pearson_binom": [], "pearson_hypergeom": []}
    pp = {
        "name": [],
        "before_pp_greatpy_size": [],
        "before_pp_local_size": [],
        "final_size": [],
        "%_of_diffrent_GO_term": [],
    }

    ann = pd.read_csv(
        "../../data/human/ontologies.csv",
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
            f"../data/tests/test_data/input/{name}", species=f"{assembly}", help=False
        )
        res_online = rpy2.robjects.r["getEnrichmentTables"](res_online)

        # local test
        # proprocessing : make a Grange frame
        df = r["read.csv"](f"../data/tests/test_data/input/{name}", sep="\t")
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
            test_file=f"../data/tests/test_data/input/{name}",
            regdom_file=f"../data/human/{assembly}/regulatory_domain.bed",
            chr_size_file=f"../data/human/{assembly}/chr_size.bed",
            annotation_file=f"../data/human/ontologies.csv",
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
