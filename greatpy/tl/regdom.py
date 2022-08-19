import pandas as pd

pd.options.display.float_format = "{:12.5e}".format


def _validate_input(association: str, max_extension: int, basal_upstream: int, basal_downstream: int) -> bool:
    """
    Checks that the inputs (association_rule, max_extension, basal_upstream, basal_downstream) are valid

    Parameters
    ----------
    association : str
        The association rule to use. Documentation aviable at https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules\n
    max_extension : int
        The maximum extension of the regulatory domain.\n
    basal_upstream : int
        The basal upstream of the regulatory domain.\n
    basal_downstream : int
        The basal downstream of the regulatory domain.\n

    Returns
    -------
    bool
        True if the inputs are valid, False otherwise.

    Examples
    --------
    >>> validate_input("two_closet")
    ...    True

    >>> validate_input("Two_Closet_Extension")
    ...    False
    """
    if association != "one_closet" and association != "two_closet" and association != "basal_plus_extention":
        print("Association rule should be OneCloset or TwoCloset Basalplusextention")
        return False
    if max_extension < 0:
        print(f"Maximum extension must be a non-negative integer: {max_extension}")
        return False
    if basal_upstream < 0:
        print(f"Basal upstream must be a non-negative integer: {basal_upstream}")
        return False
    if basal_downstream < 0:
        print(f"Basal downstream must be a non-negative integer: {basal_downstream}")
        return False
    return True


def _write_regdom(regdom: pd.DataFrame, file_name: str) -> None:
    """
    Write the regulatory regions calculated in a file given in argument

    Parameters
    ----------
    regdom : pd.DataFrame
        The regulatory regions to write.\n
    file_name : str
        The path of the file to write.\n

    Returns
    -------
    None.
    """
    f = open(file_name, "w")
    f.write("#chr\tchrstart\tchrend\tname\ttss\tstrand\n")
    for i in range(regdom.shape[0]):
        curr = regdom.iloc[i]
        chr = curr["chr"]
        start = curr["chr_start"]
        end = curr["chr_end"]
        name = curr["name"]
        tss = curr["tss"]
        strand = curr["strand"]
        f.write(f"{chr}\t{start}\t{end}\t{name}\t{tss}\t{strand}\n")
    f.close()


def _create_basal_plus_extension_regdom(
    tss: pd.DataFrame, maximumExtension: int, basalUp: int, basalDown: int, chr_size: pd.DataFrame
) -> pd.DataFrame:
    """
    Create the regulatory domains using the Basalplusextention association rule

    Parameters
    ----------
    tss : pd.DataFrame
        The TSS of the genes.\n
    maximumExtension : int
        The maximum extension of the regulatory domain.\n
    basalUp : int
        The basal upstream of the regulatory domain.\n
    basalDown : int
        The basal downstream of the regulatory domain.\n
    chr_size : pd.DataFrame
        The chromosome size.\n

    Returns
    -------
    tss : pd.DataFrame
        The regulatory domains.

    Examples
    --------
    >>> regdom = create_basal_plus_extension_regdom(
        tss_file=pd.read_csv("../../data/human/tss.bed",sep="\t",names=["chr","tss","strand"]),
        maximumExtension=100000,
        basalUp=5000,
        basalDown=1000,
        chr_sizes=pd.read_csv("../../data/human/chr_size.bed",sep="\t",names=["chr","size"])
    )
    >>> regdom.head()
    ...    |    | chr   |   chr_start |   chr_end | name      |   tss | strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     22436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       16436 |     22436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       16436 |     22436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       16436 |     28370 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       22436 |     34370 | WASH7P    | 29370 | -        |

    """
    prev = 0
    curr = 0
    next = 0
    chr_strat_end = []
    start = []
    end = []
    for i in range(tss.shape[0]):
        curr = tss.iloc[i]
        chr = curr["chr"]
        curr_chr_size = int(chr_size.loc[chr_size["chr"] == chr, "size"])
        tmp = curr["tss"]
        if curr["strand"] == "+":
            curr_chr_start = max(0, tmp - basalUp)
            curr_chr_end = min(curr_chr_size, tmp + basalDown)
        elif curr["strand"] == "-":
            curr_chr_start = max(0, tmp - basalDown)
            curr_chr_end = min(curr_chr_size, tmp + basalUp)
        elif curr["strand"] == ".":
            print("Invalid_Input : Impossible to create a basal expression regdom if you have not specify the strand")
            return False
        else:
            err = curr["strand"]
            print(f"Invalid input : strand should be '+' or '-'. Line {i} : strand = {err}")
            return False
        chr_strat_end.append([curr_chr_start, curr_chr_end])

    for i in range(tss.shape[0]):
        curr = tss.iloc[i]
        chr = curr["chr"]

        curr_chr_size = int(chr_size.loc[chr_size["chr"] == chr, "size"])
        if i != tss.shape[0] - 1:
            next = tss.iloc[i + 1]
        else:
            next = 0

        tmp_start = max(0, curr["tss"] - maximumExtension)
        basal_start = chr_strat_end[i][0]
        tmp_start = min(basal_start, tmp_start)
        if type(prev) != int and prev["chr"] == curr["chr"]:
            if prev["strand"] == "+":
                prev_end = prev["tss"] + basalDown
            else:
                prev_end = prev["tss"] + basalUp
            tmp_start = min(basal_start, max(prev_end, tmp_start))

        tmp_end = min(curr_chr_size, curr["tss"] + 1000000)
        basal_end = chr_strat_end[i][1]

        tmp_end = max(basal_end, tmp_end)
        if type(next) != int and next["chr"] == curr["chr"]:
            if next["strand"] == "+":
                nextstart = next["tss"] - basalUp
            else:
                nextstart = next["tss"] - basalDown
            tmp_end = max(basal_end, min(nextstart, tmp_end))

        prev = tss.iloc[i]
        start.append(int(tmp_start))
        end.append(int(tmp_end))
    tss["chr_start"] = start
    tss["chr_end"] = end
    return tss


def _create_two_closet_regdom(tss: pd.DataFrame, max_extension: int, chr_size: pd.DataFrame) -> pd.DataFrame:
    """
    Create the regulatory domains using the TwoCloset association rule.\n
    It is based on the basal plus extension rule but with basalUp and basalDown equals to 0.\n

    Parameters
    ----------
    tss : pd.DataFrame
        The TSS of the genes.\n
    maximumExtension : int
        The maximum extension of the regulatory domain.\n
    chr_size : pd.DataFrame
        The chromosome size.\n

    Returns
    -------
    pd.DataFrame
        The regulatory domains.

    Examples
    --------
    >>> regdom = _create_two_closet_regdom(
        tss_file=pd.read_csv("../../data/human/tss.bed",sep="\t",names=["chr","tss","strand"]),
        maximumExtension=100000,
        chr_sizes=pd.read_csv("../../data/human/chr_size.bed",sep="\t",names=["chr","size"])
    )
    >>> regdom.head()
    ...    |    | chr   |   chr_start |   chr_end | name      |   tss | strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     17436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       17436 |     17436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       17436 |     17436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       17436 |     29370 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       17436 |     30365 | WASH7P    | 29370 | -        |

    """
    return _create_basal_plus_extension_regdom(tss, max_extension, 0, 0, chr_size)


def _create_one_closet_regdom(tss: pd.DataFrame, maximum_extension: int, chr_size: pd.DataFrame) -> pd.DataFrame:
    """
    Create the regulatory domains using the OneCloset association rule

    Parameters
    ----------
    tss : pd.DataFrame
        The TSS of the genes.\n
    maximumExtension : int
        The maximum extension of the regulatory domain.\n
    basalUp : int
        The basal upstream of the regulatory domain.\n
    basalDown : int
        The basal downstream of the regulatory domain.\n
    chr_size : pd.DataFrame
        The chromosome size.\n

    Returns
    -------
    tss : pd.DataFrame
        The regulatory domains.

    Examples
    --------
    >>> regdom = create_basal_plus_extension_regdom(
        tss_file=pd.read_csv("../../data/human/tss.bed",sep="\t",names=["chr","tss","strand"]),
        maximum_extension=100000,
        chr_sizes=pd.read_csv("../../data/human/chr_size.bed",sep="\t",names=["chr","size"])
    )
    >>> regdom.head()
    ...    |    | chr   |   chr_start |   chr_end | name      |   tss | strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     17436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       17436 |     17436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       17436 |     17436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       17436 |     23403 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       23403 |     29867 | WASH7P    | 29370 | -        |

    """
    prev = 0
    curr = 0
    next = 0
    start = []
    end = []

    for i in range(tss.shape[0]):
        curr = tss.iloc[i]
        chr = curr["chr"]
        curr_chr_size = int(chr_size.loc[chr_size["chr"] == chr, "size"])
        if i < tss.shape[0] - 1:
            next = tss.iloc[i + 1]
        else:
            next = 0

        if (type(prev) == int and prev == 0) or (type(prev) != int and prev["chr"] != curr["chr"]):
            prev = 0
        if (type(next) == int and next == 0) or (type(next) != int and next["chr"] != curr["chr"]):
            next = 0

        tmp_start = max(0, curr["tss"] - maximum_extension)
        if type(prev) != int:
            middle = (curr["tss"] + prev["tss"]) // 2
            tmp_start = max(tmp_start, middle)

        tmp_end = min(curr["tss"] + maximum_extension, curr_chr_size)
        if type(next) != int:
            middle = (curr["tss"] + next["tss"]) // 2
            tmp_end = min(tmp_end, middle)

        start.append(tmp_start)
        end.append(tmp_end)
        prev = curr

    tss["chr_start"] = start
    tss["chr_end"] = end
    return tss


def create_regdom(
    tss_file: str,
    chr_sizes_file: str,
    association_rule: str,
    max_extension: int = 1000000,
    basal_upstream: int = 5000,
    basal_downstream: int = 1000,
    out_path: str or None = None,
) -> pd.DataFrame:
    """\
    Create regdoms according to the three association rules, to write the result in a file or not and to return the result as a pd.DataFrame

    Parameters
    ----------
    tss_file : str
        The path of the TSS file.\n
    chr_sizes_file : str
        The path of the chromosome size file.\n
    association_rule : str
        The association rule to use. Could be : "one_closet", "two_closet", "basal_plus_extention".\n
        Documentation aviable at https://great-help.atlassian.net/wiki/spaces/GREAT/pages/655443/Association+Rules.\n
    maximumExtension : int
        The maximum extension of the regulatory domain.\n
        Default is `100000`\n
    basalUp : int
        The basal upstream of the regulatory domain.\n
        Default is `5000`\n
    basalDown : int
        The basal downstream of the regulatory domain.\n
        Default is `1000`\n
    out_path : str or NoneType
        The path of the output file.\n
        If None, the result is only returned as a pd.DataFrame.\n
        Default is `None`\n

    Returns
    -------
    out : pd.DataFrame
        The regulatory domains.

    Examples
    --------
    >>> regdom = create_regdom(
        tss_file="../../data/human/tss.bed",
        chr_sizes_file="../../data/human/chr_size.bed",
        sep="\t",
        names=["chr","size"],
        association_rule="one_closet"
        )
    >>> regdom.head()
    ...    |    | chr   |   chr_start |   chr_end | name      |   tss | strand   |
    ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
    ...    |  0 | chr1  |           0 |     17436 | MIR6859-1 | 17436 | -        |
    ...    |  1 | chr1  |       17436 |     17436 | MIR6859-2 | 17436 | -        |
    ...    |  2 | chr1  |       17436 |     17436 | MIR6859-3 | 17436 | -        |
    ...    |  3 | chr1  |       17436 |     23403 | MIR6859-4 | 17436 | -        |
    ...    |  4 | chr1  |       23403 |     29867 | WASH7P    | 29370 | -        |

    """
    if not _validate_input(association_rule, max_extension, basal_upstream, basal_downstream):
        print("Invalid input")
        return False
    df = pd.read_csv(tss_file, sep="\t", comment="#", names=["chr", "tss", "strand", "name"])

    df = df.sort_values(["chr", "tss", "strand", "name"])

    chr_size = pd.read_csv(chr_sizes_file, sep="\t", comment="#", names=["chr", "size"])

    if association_rule == "one_closet":
        out = _create_one_closet_regdom(df, max_extension, chr_size)
    elif association_rule == "two_closet":
        out = _create_two_closet_regdom(df, max_extension, chr_size)
    elif association_rule == "basal_plus_extention":
        out = _create_basal_plus_extension_regdom(df, max_extension, basal_upstream, basal_downstream, chr_size)
    else:
        return False
    out = out.astype({"chr_start": int, "chr_end": int})
    out = out.reindex(["chr", "chr_start", "chr_end", "name", "tss", "strand"], axis=1)

    if out_path is not None:
        _write_regdom(out, out_path)
    return out
