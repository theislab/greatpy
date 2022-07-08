from types import NoneType
from anndata import AnnData
import pandas as pd
pd.options.display.float_format = '{:12.5e}'.format


def basic_tool(adata: AnnData) -> int:
    """Run a tool on the AnnData object."""
    print("Implement a tool to run on the AnnData object.")
    return 0

def validate_input(association:str,max_extension:int,basal_upstream:int,basal_downstream:int): 
    """
    This function checks that the inputs (association_rule, max_extension, basal_upstream, basal_downstream) are valid 

    Parameters
    ----------
    association : str 
        The association rule to use.
    max_extension : int
        The maximum extension of the regulatory domain.
    basal_upstream : int
        The basal upstream of the regulatory domain.
    basal_downstream : int
        The basal downstream of the regulatory domain.

    Returns
    -------
    bool
        True if the inputs are valid, False otherwise.

    Exemples 
    --------
    >>> validate_input("Two_Closet")
        True

    >>> validate_input("Two_Closet_Extension")
        False
    """
    if association != "OneCloset" and association != "TwoCloset" and association != "Basalplusextention" : 
        print("Association rule should be OneCloset or TwoCloset Basalplusextention")
        return False
    if max_extension < 0: 
        print(f"Maximum extension must be a non-negative integer: {max_extension}")
        return False
    if basal_upstream < 0 : 
        print(f"Basal upstream must be a non-negative integer: {basal_upstream}")
        return False
    if basal_downstream < 0 : 
        print(f"Basal downstream must be a non-negative integer: {basal_downstream}")
        return False
    return True

def write_Regdom(regdom:pd.DataFrame,file_name:str):
    """
    This method allows you to write the regulatory regions calculated in a file given in argument

    Parameters
    ----------
    regdom : pd.DataFrame 
        The regulatory regions to write.
    file_name : str
        The path of the file to write.

    Returns
    -------
    None.
    """
    f = open(file_name,"w")
    f.write("#chr\tChrStart\tChrEnd\tname\ttss\tstrand\n")
    for i in range(regdom.shape[0]): 
        curr = regdom.iloc[i]
        chr = curr["Chr"];start = curr["Chr_Start"];end = curr["Chr_End"];name = curr["name"];tss = curr["tss"];strand = curr["Strand"]
        f.write(f"{chr}\t{start}\t{end}\t{name}\t{tss}\t{strand}\n")
    f.close()

def create_basal_plus_extension_regdom(tss:pd.DataFrame,maximumExtension:int,basalUp:int,basalDown:int,chr_size:pd.DataFrame):
    """
    This function allows to create the regulatory domains using the Basalplusextention association rule

    Parameters
    ----------
    tss : pd.DataFrame
        The TSS of the genes.
    maximumExtension : int
        The maximum extension of the regulatory domain.
    basalUp : int
        The basal upstream of the regulatory domain.
    basalDown : int
        The basal downstream of the regulatory domain.
    chr_size : pd.DataFrame
        The chromosome size.

    Returns
    -------
    tss : pd.DataFrame
        The regulatory domains.

    Exemples 
    --------
    >>> regdom = create_basal_plus_extension_regdom(
    ...    tss_file=pd.read_csv("../../data/human/tss.bed",sep="\t",names=["Chr","tss","Strand"]),
    ...    maximumExtension=100000,
    ...    basalUp=5000,
    ...    basalDown=1000,
    ...    chr_sizes=pd.read_csv("../../data/human/chr_size.bed",sep="\t",names=["Chr","Size"])
    ... )
    >>> regdom.head()
    ...    |    | Chr   |   Chr_Start |   Chr_End | name      |   tss | Strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     22436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       16436 |     22436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       16436 |     22436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       16436 |     28370 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       22436 |     34370 | WASH7P    | 29370 | -        |
    
    """
    prev = curr = next = 0
    chr_strat_end = []
    start = []
    end = []
    for i in range(tss.shape[0]):
        curr = tss.iloc[i]
        chr = curr["Chr"]
        curr_chr_size = int(chr_size.loc[chr_size["Chr"]==chr,"Size"])
        tmp = curr["tss"]
        if curr["Strand"] == "+": 
            curr_chr_start = max(0,tmp-basalUp)
            curr_chr_end = min(curr_chr_size,tmp+basalDown)
        elif curr["Strand"] == "-": 
            curr_chr_start = max(0,tmp-basalDown)
            curr_chr_end = min(curr_chr_size,tmp+basalUp)
        elif curr["Strand"] == "." : 
            print("Invalid_Input : Impossible to create a basal expression regdom if you have not specify the strand")
            return False
        else : 
            err = curr["Strand"]
            print (f"Invalid input : strand should be '+' or '-'. Line {i} : strand = {err}")
            return False
        chr_strat_end.append([curr_chr_start,curr_chr_end])

    for i in range(tss.shape[0]):
        curr = tss.iloc[i]
        chr = curr["Chr"]
        
        curr_chr_size = int(chr_size.loc[chr_size["Chr"]==chr,"Size"])
        if i != tss.shape[0]-1 : 
            next = tss.iloc[i+1]
        else : next = 0


        tmp_start = max(0,curr["tss"]-maximumExtension)
        basal_start = chr_strat_end[i][0]
        tmp_start = min(basal_start,tmp_start)
        if type(prev) != int and prev["Chr"] == curr["Chr"]:
            if prev["Strand"] == "+": 
                prev_end = prev["tss"]+basalDown
            else : 
                prev_end = prev["tss"]+basalUp
            tmp_start = min(basal_start,max(prev_end,tmp_start))


        tmp_end = min(curr_chr_size,curr["tss"]+1000000)
        basal_end = chr_strat_end[i][1]
        
        tmp_end = max(basal_end,tmp_end)
        if type(next) != int and next["Chr"] == curr["Chr"]:
            if next["Strand"] == "+": 
                nextStart = next["tss"]-basalUp
            else : 
                nextStart = next["tss"]-basalDown
            tmp_end = max(basal_end,min(nextStart,tmp_end))
            
        prev = tss.iloc[i]
        start.append(int(tmp_start))
        end.append(int(tmp_end))
    tss["Chr_Start"] = start
    tss["Chr_End"] = end
    return tss

def create_Two_Closet_Regdom(tss:pd.DataFrame,max_extension:int,chr_size:pd.DataFrame):
    """
    This function allows to create the regulatory domains using the TwoCloset association rule. 
    It is based on the basal plus extension rule but with basalUp and basalDown equals to 0.

    Parameters
    ----------
    tss : pd.DataFrame
        The TSS of the genes.
    maximumExtension : int
        The maximum extension of the regulatory domain.
    chr_size : pd.DataFrame
        The chromosome size.

    Returns
    -------
    pd.DataFrame
        The regulatory domains.

    Exemples 
    --------
    >>> regdom = create_Two_Closet_regdom(
    ...    tss_file=pd.read_csv("../../data/human/tss.bed",sep="\t",names=["Chr","tss","Strand"]),
    ...    maximumExtension=100000,
    ...    chr_sizes=pd.read_csv("../../data/human/chr_size.bed",sep="\t",names=["Chr","Size"])
    ... )
    >>> regdom.head()
    ...    |    | Chr   |   Chr_Start |   Chr_End | name      |   tss | Strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     17436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       17436 |     17436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       17436 |     17436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       17436 |     29370 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       17436 |     30365 | WASH7P    | 29370 | -        |
    
    """
    return create_basal_plus_extension_regdom(tss,max_extension,0,0,chr_size)

def create_one_closet_regdom(tss:pd.DataFrame,maximum_extension:int,chr_size:pd.DataFrame):
    """
    This function allows to create the regulatory domains using the OneCloset association rule

    Parameters
    ----------
    tss : pd.DataFrame
        The TSS of the genes.
    maximumExtension : int
        The maximum extension of the regulatory domain.
    basalUp : int
        The basal upstream of the regulatory domain.
    basalDown : int
        The basal downstream of the regulatory domain.
    chr_size : pd.DataFrame
        The chromosome size.

    Returns
    -------
    tss : pd.DataFrame
        The regulatory domains.

    Exemples 
    --------
    >>> regdom = create_basal_plus_extension_regdom(
    ...    tss_file=pd.read_csv("../../data/human/tss.bed",sep="\t",names=["Chr","tss","Strand"]),
    ...    maximum_extension=100000,
    ...    chr_sizes=pd.read_csv("../../data/human/chr_size.bed",sep="\t",names=["Chr","Size"])
    ... )
    >>> regdom.head()
    ...    |    | Chr   |   Chr_Start |   Chr_End | name      |   tss | Strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     17436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       17436 |     17436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       17436 |     17436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       17436 |     23403 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       23403 |     29867 | WASH7P    | 29370 | -        |
    
    """
    prev = curr = next = 0
    start = []
    end = []

    for i in range(tss.shape[0]):
        curr = tss.iloc[i]
        chr = curr["Chr"]
        curr_chr_size = int(chr_size.loc[chr_size["Chr"]==chr,"Size"])
        if i < tss.shape[0]-1 : next = tss.iloc[i+1]
        else : next = 0

        if (type(prev) == int and prev == 0) or (type(prev) != int and prev["Chr"] != curr["Chr"]) : prev = 0 
        if (type(next) == int and next == 0) or (type(next) != int and next["Chr"] != curr["Chr"]): next = 0

        tmp_start = max(0,curr["tss"]-maximum_extension)
        if type(prev) != int : 
            middle = (curr["tss"]+prev["tss"])//2
            tmp_start = max(tmp_start,middle)
        
        tmp_end = min(curr["tss"]+maximum_extension,curr_chr_size)
        if type(next) != int  :
            middle = (curr["tss"]+next["tss"])//2
            tmp_end = min(tmp_end,middle)
        
        start.append(tmp_start)
        end.append(tmp_end)
        prev = curr
        
    tss["Chr_Start"] = start
    tss["Chr_End"] = end
    return tss

def create_regdom(tss_file,chr_sizes_file,association_rule,max_extension:int=1000000,basal_upstream:int=5000,basal_downstream:int=1000, out_path:str or NoneType=None): 
    """
    This function allows to create regdoms according to the three association rules, to write the result in a file or not and to return the result as a pd.DataFrame

    Parameters
    ----------
    tss_file : str
        The path of the TSS file.
    chr_sizes_file : str
        The path of the chromosome size file.
    association_rule : str
        The association rule to use. Could be : "OneCloset", "TwoCloset", "BasalPlusExtension".
    maximumExtension : int : default 1000000
        The maximum extension of the regulatory domain.
    basalUp : int : default 5000
        The basal upstream of the regulatory domain.
    basalDown : int : default 1000
        The basal downstream of the regulatory domain.
    out_path : str or NoneType : default None
        The path of the output file.

    Returns
    -------
    out : pd.DataFrame
        The regulatory domains.

    Exemples 
    --------
    >>> regdom = create_regdom(
    ...    tss_file="../../data/human/tss.bed",
    ...    chr_sizes_file="../../data/human/chr_size.bed",sep="\t",names=["Chr","Size"],
    ...    association_rule="OneCloset",
    ... )
    >>> regdom.head()
    ...    |    | Chr   |   Chr_Start |   Chr_End | name      |   tss | Strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     17436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       17436 |     17436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       17436 |     17436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       17436 |     23403 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       23403 |     29867 | WASH7P    | 29370 | -        |

    """
    if not validate_input(association_rule): 
        print("Invalid input")
        return False
    df = pd.read_csv(tss_file,sep="\t",comment="#",names=["Chr","tss","Strand","name"])
    
    df = df.sort_values(["Chr","tss","Strand","name"])

    chr_size = pd.read_csv(chr_sizes_file,sep="\t",comment="#",names=["Chr","Size"])

    if association_rule == "OneCloset" : 
        out = create_one_closet_regdom(df,max_extension,chr_size)
    elif association_rule == "TwoCloset" : 
        out = create_Two_Closet_Regdom(df,max_extension,chr_size)
    elif association_rule == "Basalplusextention" : 
        out = create_basal_plus_extension_regdom(df,max_extension,basal_upstream,basal_downstream,chr_size)
    else : 
        return False
    out = out.astype({"Chr_Start":int,"Chr_End":int})
    out = out.reindex(["Chr","Chr_Start","Chr_End","name","tss","Strand"],axis=1)

    if out_path != None : 
        write_Regdom(out,out_path) 
    return out

