import greatpy as gp

def get_nb_asso_per_region(test,regdom) : 
    """
    Function allowing from a file of genomic regions from CHIPseq 
    and a file of genomic regulatory domains to determine number of peaks 
    associated with each gene in the regulatory domain. 

    Parameters
    ----------
    test : str
        path of the file with the tests pics => columns: ["Chr","Chr_Start","Chr_End"]
    
    regdom : str
        path of the file with the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    res : dict
        dict with the number of associated genes per genomic region : key = associated gene, value = number of peaks associated with the gene 
        
    Exemples 
    --------
    test = pd.DataFrame({
    ...    "Chr":["chr1"],
    ...    "Chr_Start":[1052028],
    ...    "Chr_End": [1052049]})

    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']})

    >>> get_association(test,regdom)        
        {'RNF223':2}
    
    """
    res = {}
    test,regdom,_,_ = gp.tl.GREAT.loader(test,regdom,None,None) 
    for i in range(test.shape[0]) :
        currTest = test.iloc[i]
        regdom_curr_test = regdom.loc[(regdom["Chr"] == currTest["Chr"])].sort_values("Chr_Start")
        regdom_curr_test = regdom_curr_test.loc[
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"])) | # regdom overlap totally test 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"])) | # test overlap totally regdom 
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_Start"])) | # regdom overlap not totally test on left side 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"]) & (regdom_curr_test["Chr_Start"] <= currTest["Chr_End"])) # regdom overlap not totally test on right side 
            ] 
        res[i] = regdom_curr_test.shape[0]
    return res

def get_dist_to_tss(test,regdom) : 
    """
    Function allowing from a file of genomic regions from CHIPseq 
    and a file of genomic regulatory domains to determine the distance from peaks 
    to the transcription start site of the associated gene

    Parameters
    ----------
    test : str
        path of the file with the tests pics => columns: ["Chr","Chr_Start","Chr_End"]
    
    regdom : str
        path of the file with the regulatory domains => columns: ["Chr"	"Chr_Start"	"Chr_End"	"Name"	"tss"	"strand"].

    Returns
    -------
    res : dict
        dict with the distance from tss to the associated genes : key = number of the input, value = distance from peaks to tss of associated genes 
        
    Exemples 
    --------
    test = pd.DataFrame({
    ...    "Chr":["chr1"],
    ...    "Chr_Start":[1052028],
    ...    "Chr_End": [1052049]})

    regdom = pd.DataFrame({
    ...    "Chr":["chr1","chr1"],
    ...    "Chr_Start":[1034992,1079306],
    ...    "Chr_End": [1115089,1132016],
    ...    "Name":["RNF223","C1orf159"],
    ...    "tss":[1074306,1116089],
    ...    "strand":['-','-']})

    >>> get_association(test,regdom)        
        {'RNF223':[-22278]}
    
    """
    res = {}
    test,regdom,_,_ = gp.tl.GREAT.loader(test,regdom,None,None) 
    for i in range(test.shape[0]) :
        currTest = test.iloc[i]
        mean_pos_test = (currTest["Chr_End"] + currTest["Chr_Start"])/2
        regdom_curr_test = regdom.loc[(regdom["Chr"] == currTest["Chr"])].sort_values("Chr_Start")
        regdom_curr_test = regdom_curr_test.loc[
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"])) | # regdom overlap totally test 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"])) | # test overlap totally regdom 
            ((regdom_curr_test["Chr_Start"] <= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] <= currTest["Chr_End"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_Start"])) | # regdom overlap not totally test on left side 
            ((regdom_curr_test["Chr_Start"] >= currTest["Chr_Start"]) & (regdom_curr_test["Chr_End"] >= currTest["Chr_End"]) & (regdom_curr_test["Chr_Start"] <= currTest["Chr_End"])) # regdom overlap not totally test on right side 
            ] 
        res[i] = []
        for j in range(regdom_curr_test.shape[0]) : 
            res[i].append(regdom_curr_test.iloc[j]["tss"] - mean_pos_test)
    return res