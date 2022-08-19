import pandas as pd

import greatpy as great


def get_nb_asso_per_region(test: str or pd.DataFrame, regdom: str or pd.DataFrame) -> dict:
    """
    Determine number of peaks associated with each gene in the regulatory domain.

    Parameters
    ----------
    test : str or pd.DataFrame
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]\n
    regdom : str or pd.DataFrame
        path of the file with the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"name"	"tss"	"strand"].\n

    Returns
    -------
    res : dict
        dict with the number of associated genes per genomic region : \n
        - key = associated gene\n
        - value = number of peaks associated with the gene

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
    test, regdom, _, _ = great.tl.loader(test, regdom, None, None)
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
        path of the file with the tests pics => columns: ["chr","chr_start","chr_end"]\n
    regdom : str or pd.DataFrame
        path of the file with the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"name"	"tss"	"strand"].\n

    Returns
    -------
    res : dict
        dict with the distance from tss to the associated genes : \n
        - key = number of the input \n
        - value = distance from peaks to tss of associated genes\n

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
    test, regdom, _, _ = great.tl.loader(test, regdom, None, None)
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
