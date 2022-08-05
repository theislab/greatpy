from math import exp, fabs, lgamma, log

import dask.dataframe as dd
import numpy as np
import pandas as pd
from scipy.special import comb
from statsmodels.stats.multitest import multipletests
import bindome as bd 
pd.options.display.float_format = "{:12.5e}".format


class GREAT:
    def loader(
        test_data: None or str or pd.DataFrame,
        regdom_file: None or str or pd.DataFrame,
        chr_size_file: None or str or pd.DataFrame,
        annotation_file: None or str or pd.DataFrame,
    ):
        """
        Used to load all datasets needed for the enrichment calculation

        Parameters
        ----------
        test_data : None or str or pd.DataFrame
            Genomic set of peaks to be tested
        regdom_file : None or str or pd.DataFrame
            Regulatory domain of all genes in the genome
        chr_size_file : None or str or pd.DataFrame
            Table with the size of each chromosome
        annotation_file : None or str or pd.DataFrame
            Table with the annotation of each gene in the genome

        Returns
        -------
        test_data : pd.DataFrame
            Genomic set of peaks to be tested in the good format
        regdom : pd.DataFrame
            Regulatory domain of all genes in the genome in the good format
        size : pd.DataFrame
            Table with the size of each chromosome in the good format
        ann : pd.DataFrame
            Table with the annotation of each gene in the genome in the good format

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../../data/tests/test_data/input/02_srf_hg38.bed",
            "../../data/human/hg38/regulatory_domain.bed",
            "../../data/human/hg38/chr_size.bed",
            "../data/human/ontologies.csv"
            )

        >>> test.head()
        ...    |    | Chr   |   Chr_Start |   Chr_End |
        ...    |---:|:------|------------:|----------:|
        ...    |  0 | chr1  |     1052028 |   1052049 |
        ...    |  1 | chr1  |     1065512 |   1065533 |
        ...    |  2 | chr1  |     1067375 |   1067397 |
        ...    |  3 | chr1  |     1068083 |   1068119 |
        ...    |  4 | chr1  |    10520283 |  10520490 |

        >>> regdom.head()
        ...    |    | Chr   |   Chr_Start |   Chr_End | Name      |   tss | Strand   |
        ...    |---:|:------|------------:|----------:|:----------|------:|:---------|
        ...    |  0 | chr1  |           0 |     22436 | MIR6859-1 | 17436 | -        |
        ...    |  1 | chr1  |       16436 |     22436 | MIR6859-2 | 17436 | -        |
        ...    |  2 | chr1  |       16436 |     22436 | MIR6859-3 | 17436 | -        |
        ...    |  3 | chr1  |       16436 |     28370 | MIR6859-4 | 17436 | -        |
        ...    |  4 | chr1  |       22436 |     34370 | WASH7P    | 29370 | -        |

        >>> size.head()
        ...    |    | Chrom   |      Size |
        ...    |---:|:--------|----------:|
        ...    |  0 | chr1    | 248956422 |
        ...    |  1 | chr2    | 242193529 |
        ...    |  2 | chr3    | 198295559 |
        ...    |  3 | chr4    | 190214555 |
        ...    |  4 | chr5    | 181538259 |

        >>> ann.head()
        ...    |    | id         | name                                                   | symbol        |
        ...    |---:|:-----------|:-------------------------------------------------------|:--------------|
        ...    |  0 | GO:0003924 | GTPase activity                                        | DNAJC25-GNG10 |
        ...    |  1 | GO:0007186 | G protein-coupled receptor signaling pathway           | DNAJC25-GNG10 |
        ...    |  2 | GO:0003723 | RNA binding                                            | NUDT4B        |
        ...    |  3 | GO:0005829 | cytosol                                                | NUDT4B        |
        ...    |  4 | GO:0008486 | diphosphoinositol-polyphosphate diphosphatase activity | NUDT4B        |

        """

        if type(regdom_file) == str:
            regdom = pd.read_csv(
                regdom_file,
                sep="\t",
                comment="#",
                names=["chr", "chr_start", "chr_end", "Name", "tss", "Strand"],
                dtype={
                    "chr": "object",
                    "chr_start": "int64",
                    "chr_end": "int64",
                    "Name": "object",
                    "tss": "int64",
                    "Strand": "object",
                },
            )
        elif type(regdom_file) == pd.DataFrame:
            regdom = regdom_file.iloc[:, :6]
            colname = list(regdom.columns)
            try:
                regdom = regdom.rename(
                    columns={
                        colname[0]: "chr",
                        colname[1]: "chr_start",
                        colname[2]: "chr_end",
                        colname[3]: "Name",
                        colname[4]: "tss",
                        colname[5]: "Strand",
                    }
                )
            except:
                print("Error in the format of the regdom file")
                print("The regdom file must have the following columns : Chr, Chr_Start, Chr_End, Name, tss, Strand")
                return False
        else:
            regdom = regdom_file

        if type(test_data) == str:
            test_data = pd.read_csv(
                test_data,
                sep="\t",
                comment="#",
                usecols=[0, 1, 2],
                names=["chr", "chr_start", "chr_end"],
                dtype={"chr": "object", "chr_start": "int64", "chr_end": "int64"},
            )
        elif type(test_data) == pd.DataFrame:
            test_data = test_data.iloc[:, :3]
            colname = list(test_data.columns)
            try:
                test_data = test_data.rename(
                    columns={colname[0]: "chr", colname[1]: "chr_start", colname[2]: "chr_end"}
                )
            except:
                print("Error in test dataframe, please check your input")
                print("Columns should be : chr...(type object), start(type int), end(type int)")
                return False
        else:
            pass

        if type(chr_size_file) == str:
            size = pd.read_csv(
                chr_size_file,
                sep="\t",
                comment="#",
                names=["Chrom", "Size"],
                dtype={"Chrom": "object", "Size": "int64"},
            )
        elif type(chr_size_file) == pd.DataFrame:
            size = chr_size_file.iloc[:, :2]
            colname = list(size.columns)
            try:
                size = size.rename(columns={colname[0]: "Chrom", colname[1]: "Size"})
            except:
                print("Error in the format of the chr_size file")
                print("The chr_size file must have the following columns : Chrom, Size")
                return False
        else:
            size = chr_size_file

        if type(annotation_file) == str:
            dask_df = dd.read_csv(
                annotation_file,
                sep=";",
                comment="#",
                dtype={
                    "ensembl": "object",
                    "id": "object",
                    "name": "object",
                    "ontology.group": "object",
                    "gene.name": "object",
                    "symbol": "object",
                },
                usecols=["id", "name", "symbol"],
                low_memory=False,
            )
            ann = dask_df.compute()
            ann = ann[ann["id"].str.match("^GO.*") == True]
        elif type(annotation_file) == pd.DataFrame:
            ann = annotation_file.iloc[:, :4]
            colname = list(ann.columns)
            try:
                ann = ann.rename(columns={colname[0]: "id", colname[1]: "name", colname[3]: "symbol"})
            except:
                print("Error in the format of the annotation file")
                print("The annotation file must have the following columns : id, name, symbol")
                return False
        else:
            ann = annotation_file

        return test_data, regdom, size, ann

    def __enrichment_binom_and_hypergeom(
        test: pd.DataFrame, regdom: pd.DataFrame, size: pd.DataFrame, ann: pd.DataFrame, asso: list
    ) -> pd.DataFrame:
        """
        Used to compute the enrichment of the test data using the binomial test and the hypergeometric test.

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame
            Regulatory domain of all genes in the genome
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the binomial test and the the hypergeometric test

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../../data/tests/test_data/input/02_srf_hg38.bed",
            "../../data/human/hg38/regulatory_domain.bed",
            "../../data/human/hg38/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = GREAT.____enrichment_binom_and_hypergeom(
            test = test,
            regdom = regdom,
            size = size,
            ann = ann,
            asso = get_association(test,regdom)
            )

        >>> enrichment.head()
        ...    |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |
        ...    |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|
        ...    | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  |
        ...    | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  |
        ...    | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  |
        ...    | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 |
        ...    | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  |
        """
        # Init Great
        res = {}
        hit = {}

        # init Hypergeom
        hypergeom_gene_set = len(asso)  # get the number of genes in the test gene set.
        hypergeom_total_number_gene = regdom.shape[0]  # get the number of genes in the genome.

        # Init binom
        n_binom = test.shape[0]  # get the number of genomic region in the test set
        total_nu = size["Size"].sum()  # get the total number of nucleotides in the genome

        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[
            regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))
        ]  # reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association
        len_on_chr = len_regdom(regdom)  # get the length of each regulatory domain

        # Compute for all associating gene and for each GO id associated with the gene the probability.
        for name in asso:
            ann_name_gene = ann[ann["symbol"].isin([name])]
            id = ann_name_gene["id"]
            ann_reduce = ann[ann["id"].isin(list(id))]
            tmp = []
            for i in list(id.unique()):
                gene_imply = ann_reduce[ann_reduce["id"].isin([i])]
                K_hypergeom = gene_imply.shape[0]  # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[
                    0
                ]  # get the number of genes in the test gene set with annotation

                if i not in list(hit.keys()):
                    hit[i] = number_of_hits(
                        test, curr_regdom
                    )  # get the number of test genomic regions in the regulatory domain of a gene with annotation
                k_binom = hit[i]
                nb_binom = sum(
                    len_on_chr[i] for i in curr_regdom["Name"]
                )  # get the portion of the genome in the regulatory domain of a gene with annotation
                tmp.append((k_binom, nb_binom, i, gene_imply.iloc[0]["name"], K_hypergeom, k_hypergeom))

            res.update(
                {
                    elem[2]: [
                        elem[3],
                        get_binom_pval(n_binom, elem[0], elem[1] / total_nu),
                        elem[0] / (elem[1] / total_nu),  # binom enrichment
                        hypergeom_cdf(hypergeom_total_number_gene, elem[4], hypergeom_gene_set, elem[5]),
                        (elem[5] * hypergeom_total_number_gene)
                        / (hypergeom_gene_set * elem[4]),  # Hypergeom enrichment
                        elem[0],
                        elem[0] / elem[4],
                    ]
                    for elem in tmp
                }
            )

        return (
            pd.DataFrame(res)
            .transpose()
            .rename(
                columns={
                    0: "go_term",
                    1: "binom_p_value",
                    2: "binom_fold_enrichment",
                    3: "hypergeom_p_value",
                    4: "hypergeometric_fold_enrichment",
                    5: "intersection_size",
                    6: "recall",
                }
            )
            .replace(0, np.nan)
            .dropna()
            .sort_values(by="binom_p_value")
        )

    def __enrichment_binom(
        test: pd.DataFrame, regdom: pd.DataFrame, size: pd.DataFrame, ann: pd.DataFrame, asso: list
    ) -> pd.DataFrame:
        """
        Used to compute the enrichment of the test data using the binomial test.

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame
            Regulatory domain of all genes in the genome
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the binomial test

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../../data/tests/test_data/input/02_srf_hg38.bed",
            "../../data/human/hg38/regulatory_domain.bed",
            "../../data/human/hg38/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = GREAT.____enrichment_binom(
            test = test,
            regdom = regdom,
            size = size,
            ann = ann,
            asso = get_association(test,regdom)
            )
        >>> enrichment.head()
        ...    |            | go_term                                                   |   binom_p_value |   binom_fold_enrichment |   intersection_size |   recall |
        ...    |:-----------|:----------------------------------------------------------|----------------:|------------------------:|--------------------:|---------:|
        ...    | GO:0072749 | cellular response to cytochalasin B                       |     2.21968e-12 |                227251   |                   5 |  5       |
        ...    | GO:0051623 | positive regulation of norepinephrine uptake              |     2.21968e-12 |                227251   |                   5 |  5       |
        ...    | GO:0098973 | structural constituent of postsynaptic actin cytoskeleton |     2.1174e-10  |                 91052.6 |                   5 |  1.25    |
        ...    | GO:0097433 | dense body                                                |     6.40085e-10 |                 16061.8 |                   8 |  1.33333 |
        ...    | GO:0032796 | uropod organization                                       |     2.6988e-09  |                 54544.9 |                   5 |  2.5     |
        """
        # Init Great
        res = {}
        hit = {}

        # Init binom
        n_binom = test.shape[0]  # get the number of genomic region in the test set
        total_nu = size["Size"].sum()  # get the total number of nucleotides in the genome

        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[
            regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))
        ]  # reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association
        len_on_chr = len_regdom(regdom)  # get the length of each regulatory domain

        # Compute for all associating gene and for each GO id associated with the gene the probability.
        for name in asso:
            ann_name_gene = ann[ann["symbol"].isin([name])]
            id = ann_name_gene["id"]
            ann_reduce = ann[ann["id"].isin(list(id))]
            tmp = []
            for i in list(id.unique()):
                gene_imply = ann_reduce[ann_reduce["id"].isin([i])]
                K = gene_imply.shape[0]
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]

                if i not in list(hit.keys()):
                    hit[i] = number_of_hits(
                        test, curr_regdom
                    )  # get the number of test genomic regions in the regulatory domain of a gene with annotation
                k_binom = hit[i]
                nb_binom = sum(
                    len_on_chr[i] for i in curr_regdom["Name"]
                )  # get the portion of the genome in the regulatory domain of a gene with annotation
                tmp.append(k_binom, nb_binom, i, gene_imply.iloc[0]["name"], K)

            res.update(
                {
                    elem[2]: [
                        elem[3],
                        get_binom_pval(n_binom, elem[0], elem[1] / total_nu),
                        elem[0] / (elem[1] / total_nu),  # binom enrichment
                        elem[0],
                        elem[0] / elem[4],
                    ]
                    for elem in tmp
                }
            )

        return (
            pd.DataFrame(res)
            .transpose()
            .rename(
                columns={
                    0: "go_term",
                    1: "binom_p_value",
                    2: "binom_fold_enrichment",
                    3: "intersection_size",
                    4: "recall",
                }
            )
            .sort_values(by="binom_p_value")
        )

    def __enrichment_hypergeom(test: pd.DataFrame, regdom: pd.DataFrame, ann: pd.DataFrame, asso: list) -> pd.DataFrame:
        """
        Used to compute the enrichment of the test data using the hypergeometric test.

        Parameters
        ----------
        test : pd.DataFrame
            Genomic set of peaks to be tested
        regdom : pd.DataFrame
            Regulatory domain of all genes in the genome
        chr_size :  pd.DataFrame
            Table with the size of each chromosome
        annotation : pd.DataFrame
            Table with the annotation of each gene in the genome
        asso : list
            List of the association between gene from regdom and peaks from test

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the hypergeometric test

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = GREAT.____enrichment_hypergeom(
            test = test,
            regdom = regdom,
            ann = ann,
            asso = get_association(test,regdom)
            )
        >>> enrichment.head()
        ...    |            | go_term                                                                                    |   hypergeom_p_value |   hypergeometric_fold_enrichment |   intersection_size |    recall |
        ...    |:-----------|:-------------------------------------------------------------------------------------------|--------------------:|---------------------------------:|--------------------:|----------:|
        ...    | GO:0015629 | actin cytoskeleton                                                                         |         2.25347e-06 |                          2.73071 |                  27 | 0.116883  |
        ...    | GO:1903979 | negative regulation of microglial cell activation                                          |         0.000302551 |                         17.522   |                   3 | 0.75      |
        ...    | GO:1902626 | assembly of large subunit precursor of preribosome                                         |         0.000302551 |                         17.522   |                   3 | 0.75      |
        ...    | GO:0001077 | proximal promoter DNA-binding transcription activator activity, RNA polymerase II-specific |         0.000504006 |                          1.94689 |                  29 | 0.0833333 |
        ...    | GO:0000977 | RNA polymerase II regulatory region sequence-specific DNA binding                          |         0.000511704 |                          2.03154 |                  26 | 0.0869565 |
        """
        # Init Great
        res = {}

        # Init hypergeom
        hypergeom_total_number_gene = regdom.shape[0]  # get the number of genes in the genome
        hypergeom_gene_set = len(asso)  # get the number of genes in the test gene set.

        ann_red = ann[ann["symbol"].isin(asso)]
        regdom = regdom[
            regdom["Name"].isin(list(ann[ann["id"].isin(list(ann_red["id"]))]["symbol"]))
        ]  # reduction of the regdom file by selecting only the genes whose GO ID is owned by a gene of the association

        # Compute for all associating gene and for each GO id associated with the gene the probability.
        for name in asso:
            ann_name_gene = ann[ann["symbol"] == name]
            id = ann_name_gene["id"]
            ann_reduce = ann[ann["id"].isin(list(id))]
            tmp = []
            for i in list(id.unique()):
                gene_imply = ann_reduce[ann_reduce["id"] == i]
                K_hypergeom = gene_imply.shape[0]  # get be the number of genes in the genome with annotation
                curr_regdom = regdom.loc[regdom["Name"].isin(list(gene_imply["symbol"]))]
                k_hypergeom = curr_regdom.loc[curr_regdom["Name"].isin(asso)].shape[
                    0
                ]  # get the number of genes in the test gene set with annotation
                tmp.append((i, gene_imply.iloc[0]["name"], K_hypergeom, k_hypergeom))

            res.update(
                {
                    elem[0]: [
                        elem[1],
                        hypergeom_cdf(hypergeom_total_number_gene, elem[2], hypergeom_gene_set, elem[3]),
                        (elem[3] * hypergeom_total_number_gene)
                        / (hypergeom_gene_set * elem[2]),  # hypergeom enrichment
                        elem[3],
                        elem[3] / elem[2],
                    ]
                    for elem in tmp
                }
            )

        return (
            pd.DataFrame(res)
            .transpose()
            .rename(
                columns={
                    0: "go_term",
                    1: "hypergeom_p_value",
                    2: "hypergeometric_fold_enrichment",
                    3: "intersection_size",
                    4: "recall",
                }
            )
            .replace(0, np.nan)
            .dropna()
            .sort_values(by="hypergeom_p_value")
        )

    def enrichment(
        test_file: str or pd.DataFrame,
        regdom_file: str or pd.DataFrame,
        chr_size_file: str or pd.DataFrame,
        annotation_file: str or pd.DataFrame,
        binom=True,
        hypergeom=True,
    ) -> pd.DataFrame:
        """
        Wrapper of the 3 private methods:
        * GREAT.__enrichment_binom_and_hypergeom
        * GREAT.__enrichment_binom
        * GREAT.__enrichment_hypergeom

        Parameters
        ----------
        test_file : str or pd.DataFrame
            Genomic set of peaks to be tested
        regdom_file : str or pd.DataFrame
            Regulatory domain of all genes in the genome
        chr_size_file : str or pd.DataFrame
            Table with the size of each chromosome
        annotation_file : str or pd.DataFrame
            Table with the annotation of each gene in the genome
        binom : bool (default True)
            If True, the binomial test is used.
        hypergeom : bool (default True)
            If True, the hypergeometric test is used.

        Returns
        -------
        pd.DataFrame
            dataframe contains for every GO ID associate with a every associated gene the p-value for the hypergeometric test

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = GREAT.enrichment(
            test = test,
            regdom = regdom,
            chr_size_file = size,
            ann = ann,
            binom=True,
            hypergeom=True
            )
        >>> enrichment.head()
        ...    |            | go_term                                                   |   binom_p_value |   binom_fold_enrichment |   hypergeom_p_value |   hypergeometric_fold_enrichment |   intersection_size |   recall |
        ...    |:-----------|:----------------------------------------------------------|----------------:|------------------------:|--------------------:|---------------------------------:|--------------------:|---------:|
        ...    | GO:0072749 | cellular response to cytochalasin B                       |     2.21968e-12 |                227251   |          0.0428032  |                         23.3627  |                   5 |  5       |
        ...    | GO:0051623 | positive regulation of norepinephrine uptake              |     2.21968e-12 |                227251   |          0.0428032  |                         23.3627  |                   5 |  5       |
        ...    | GO:0098973 | structural constituent of postsynaptic actin cytoskeleton |     2.1174e-10  |                 91052.6 |          0.160543   |                          5.84068 |                   5 |  1.25    |
        ...    | GO:0097433 | dense body                                                |     6.40085e-10 |                 16061.8 |          0.00141783 |                         11.6814  |                   8 |  1.33333 |
        ...    | GO:0032796 | uropod organization                                       |     2.6988e-09  |                 54544.9 |          0.00182991 |                         23.3627  |                   5 |  2.5     |

        >>> enrichment = GREAT.enrichment(
            test = test,
            regdom = regdom,
            ann = ann,
            asso = get_association(test,regdom),
            binom=True,
            hypergeom=False
            )
        ...    |            | go_term                                                   |   binom_p_value |   binom_fold_enrichment |   intersection_size |   recall |
        ...    |:-----------|:----------------------------------------------------------|----------------:|------------------------:|--------------------:|---------:|
        ...    | GO:0072749 | cellular response to cytochalasin B                       |     2.21968e-12 |                227251   |                   5 |  5       |
        ...    | GO:0051623 | positive regulation of norepinephrine uptake              |     2.21968e-12 |                227251   |                   5 |  5       |
        ...    | GO:0098973 | structural constituent of postsynaptic actin cytoskeleton |     2.1174e-10  |                 91052.6 |                   5 |  1.25    |
        ...    | GO:0097433 | dense body                                                |     6.40085e-10 |                 16061.8 |                   8 |  1.33333 |
        ...    | GO:0032796 | uropod organization                                       |     2.6988e-09  |                 54544.9 |                   5 |  2.5     |

        >>> enrichment = GREAT.enrichment(
            test = test,
            regdom = regdom,
            ann = ann,
            asso = get_association(test,regdom),
            binom=False,
            hypergeom=True
            )
        >>> enrichment.head()
        ...    |            | go_term                                                                                    |   hypergeom_p_value |   hypergeometric_fold_enrichment |   intersection_size |    recall |
        ...    |:-----------|:-------------------------------------------------------------------------------------------|--------------------:|---------------------------------:|--------------------:|----------:|
        ...    | GO:0015629 | actin cytoskeleton                                                                         |         2.25347e-06 |                          2.73071 |                  27 | 0.116883  |
        ...    | GO:1903979 | negative regulation of microglial cell activation                                          |         0.000302551 |                         17.522   |                   3 | 0.75      |
        ...    | GO:1902626 | assembly of large subunit precursor of preribosome                                         |         0.000302551 |                         17.522   |                   3 | 0.75      |
        ...    | GO:0001077 | proximal promoter DNA-binding transcription activator activity, RNA polymerase II-specific |         0.000504006 |                          1.94689 |                  29 | 0.0833333 |
        ...    | GO:0000977 | RNA polymerase II regulatory region sequence-specific DNA binding                          |         0.000511704 |                          2.03154 |                  26 | 0.0869565 |

        """
        if not binom and not hypergeom:
            return False

        test, regdom, size, ann = GREAT.loader(test_file, regdom_file, chr_size_file, annotation_file)
        asso = get_association(
            test, regdom
        )  # get the name of the regulatory domain associated to each genomic region in the test set

        if binom and hypergeom:
            return GREAT.__enrichment_binom_and_hypergeom(test, regdom, size, ann, asso)

        elif binom:
            return GREAT.__enrichment_binom(test, regdom, size, ann, asso)

        else:
            return GREAT.__enrichment_hypergeom(test, regdom, ann, asso)
    
    def enrichment_multiple(
        tests: list ,
        regdom_file: str or pd.DataFrame,
        chr_size_file: str or pd.DataFrame,
        annotation_file: str or pd.DataFrame,
        annpath: str = "../../annotation/",
        binom=True,
        hypergeom=True,
    ) -> dict:
        """
        Compute the enrichment of each GO term on multiple tests sets using bindome. 

        Parameters
        ----------
        tests : list
            List of complete name of data to compute with bindome
        regdom_file : str or pd.DataFrame
            Regulatory domain of all genes in the genome
        chr_size_file : str or pd.DataFrame
            Table with the size of each chromosome
        annotation_file : str or pd.DataFrame
            Table with the annotation of each gene in the genome
        binom : bool (default True)
            If True, the binomial test is used.
        hypergeom : bool (default True)
            If True, the hypergeometric test is used.

        Returns
        -------
        list of pd.DataFrame
            List of dataframe with the enrichment of each test

        Examples
        --------
        >>> tests = ["MAX:K-562,WA01,HeLa-S3", "BACH1:A-549,GM12878"]
        >>> enrichment = GREAT.enrichment_multiple(
                tests = tests,
                regdom_file="../data/human/hg38/regulatory_domain.bed",
                chr_size_file="../data/human/hg38/chr_size.bed",
                annotation_file="../data/human/ontologies.csv",
                binom=True,
                hypergeom=True,
            )
        >>> enrichment
        ...    {'MAX': pd.DataFrame, 
        ...    'BACH1': pd.DataFrame}
 
        """
        if not binom and not hypergeom:
            return False
        
        _, regdom, size, ann = GREAT.loader(None, regdom_file, chr_size_file, annotation_file)
        
        bd.bindome.constants.ANNOTATIONS_DIRECTORY = annpath

        res = {}
        
        for name in tests : 
            name_TF = name.split(":")[0]
            tmp_df = bd.bindome.datasets.REMAP2020.get_remap_peaks(name_TF)
            tmp = tmp_df[tmp_df[3]==name].iloc[:,0:3]
            tmp = tmp.rename(columns={"chr":'chr','start':"chr_start",'end':"chr_end"})

            asso = get_association(tmp, regdom )  # get the name of the regulatory domain associated to each genomic region in the test set

            if binom and hypergeom:
                enrichment = GREAT.__enrichment_binom_and_hypergeom(tmp, regdom, size, ann, asso)

            elif binom:
                enrichment = GREAT.__enrichment_binom(tmp, regdom, size, ann, asso)

            else:
                enrichment = GREAT.__enrichment_hypergeom(tmp, regdom, ann, asso)
            
            res[name_TF] = enrichment

        return res

    def set_bonferroni(self, alpha: float = 0.05) -> pd.DataFrame:
        """
        Create new columns in the dataframe with the Bonferroni correction

        Parameters
        ----------
        alpha : float
            alpha value for the Bonferroni correction

        Returns
        -------
        pd.DataFrame
            dataframe new columns with the Bonferroni correction for each p-value

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = great.tl.GREAT.enrichment(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv",
            binom=True,
            hypergeom=True
            )
        >>> bonferroni = GREAT.set_bonferroni(enrichment,alpha=0.05)
        >>> bonferroni.head()
        ...    |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |   binom_bonferroni |   hypergeom_bonferroni |
        ...    |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|-------------------:|-----------------------:|
        ...    | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  |        3.0754e-10  |                      1 |
        ...    | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  |        2.87384e-07 |                      1 |
        ...    | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  |        2.87384e-07 |                      1 |
        ...    | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 |        3.74892e-07 |                      1 |
        ...    | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  |        7.68931e-07 |                      1 |

        """
        for col in self.columns:
            if col in ["binom_p_value", "hypergeom_p_value"]:
                col_split = col.split("_")
                self[f"{col_split[0]}_bonferroni"] = multipletests(self[col], alpha=alpha, method="bonferroni")[1]
        return self

    def set_fdr(self, alpha: float = 0.05) -> pd.DataFrame:
        """
        Create new columns in the dataframe with the fdr correction

        Parameters
        ----------
        alpha : float
            alpha value for the fdr correction

        Returns
        -------
        pd.DataFrame
            dataframe new columns with the fdr correction for each p-value

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = great.tl.GREAT.enrichment(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv",
            binom=True,
            hypergeom=True
            )
        >>> fdr = GREAT.set_fdr(enrichment,alpha=0.05)
        >>> fdr.head()
        ...    |            | go_term                                                          |   binom_p_value |   hypergeom_p_value |   binom_fdr |   hypergeom_fdr |
        ...    |:-----------|:-----------------------------------------------------------------|----------------:|--------------------:|------------:|----------------:|
        ...    | GO:0045887 | positive regulation of synaptic growth at neuromuscular junction |     5.17744e-13 |          0.0029275  | 3.0754e-10  |       0.0913909 |
        ...    | GO:0044721 | protein import into peroxisome matrix, substrate release         |     4.83812e-10 |          0.0029275  | 9.3723e-08  |       0.0913909 |
        ...    | GO:0036250 | peroxisome transport along microtubule                           |     4.83812e-10 |          0.0029275  | 9.3723e-08  |       0.0913909 |
        ...    | GO:0016561 | protein import into peroxisome matrix, translocation             |     6.31131e-10 |          0.00584656 | 9.3723e-08  |       0.0913909 |
        ...    | GO:0047485 | protein N-terminus binding                                       |     1.2945e-09  |          0.0050377  | 1.53786e-07 |       0.0913909 |

        """
        for col in self.columns:
            if col in ["binom_p_value", "hypergeom_p_value"]:
                col_split = col.split("_")
                self[f"{col_split[0]}_fdr"] = multipletests(self[col], alpha=alpha, method="fdr_bh")[1]
        return self

    def set_threshold(self, colname: str, alpha: int = 0.05) -> pd.DataFrame:
        """
        Delete rows according to the p-value of the column taken as argument. By default the alpha value is 0.05

        Parameters
        ----------
        alpha : float
            alpha value for the fdr correction

        Returns
        -------
        pd.DataFrame
            dataframe with the rows deleted according to the p-value threshold

        Examples
        --------
        >>> test,regdom,size,ann = GREAT.loader(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv"
            )
        >>> enrichment = great.tl.GREAT.enrichment(
            "../data/tests/test_data/input/03_srf_hg19.bed",
            "../data/human/hg19/regulatory_domain.bed",
            "../data/human/hg19/chr_size.bed",
            "../data/human/ontologies.csv",
            binom=True,
            hypergeom=True
            )
        >>> enrichment.shape[0]
        ...    594

        >>> significant = GREAT.set_threshold(enrichment,colname="binom_p_value",alpha=0.05)
        >>> significant.shape[0]
        ...    310

        """
        if colname in self.columns:
            self = self.loc[self[colname] <= alpha]
        return self


######################################################################
################# Utils function used by GREAT class #################
######################################################################
def get_association(test, regdom) -> list:
    """
    From a file of genomic regions from CHIPseq
    and a file of genomic regulatory domains determine the names
    of genes associated with at least one genomic region

    Parameters
    ----------
    test : pd.dataFrame
        df of the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : pd.dataFrame
        df of the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"Name"	"tss"	"strand"].

    Returns
    -------
    res : list
        list of gene associated with at least with one test peak

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
    ...    ['RNF223']

    """
    res = []
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
        res = res + list(regdom_curr_test["Name"])
    return list(dict.fromkeys(res))


def len_regdom(regdom: pd.DataFrame) -> dict:
    """
    Calculate for each gene name from regdom the
    size of the regulatory region for this gene in the genome

    Parameters
    ----------
    regdom : pd.dataFrame
        df of the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"Name"	"tss"	"strand"].

    Returns
    -------
    dict
        dictionary in which each key corresponds to a gene name
        from regdom and the value is the size of the regulatory
        region for that gene

    Examples
    --------
    >>> regdom = pd.DataFrame(
        {
            "chr":["chr1","chr1"],
            "chr_start":[1034992,1079306],
            "chr_end": [1115089,1132016],
            "Name":["RNF223","C1orf159"],
            "tss":[1074306,1116089],
            "strand":['-','-']
            })
    >>> len_regdom(regdom)
    ...    {'RNF223': 80097, 'C1orf159': 52710}

    """
    test = regdom["chr_end"] - regdom["chr_start"]
    return pd.DataFrame({"len": list(test)}, index=regdom["Name"]).to_dict()["len"]


def number_of_hits(test, regdom) -> int:
    """
    Calculate the number of hits from several
    genomic regions and the file describing the regulatory regions

    Parameters
    ----------
    test : pd.dataFrame
        df of the tests pics => columns: ["chr","chr_start","chr_end"]

    regdom : pd.dataFrame
        df of the regulatory domains => columns: ["chr"	"chr_start"	"chr_end"	"Name"	"tss"	"strand"].

    Returns
    -------
    nb : int
        number of hit

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
    >>> number_of_hits(test,regdom)
    ...    1

    """
    nb = 0
    regdom = regdom[["chr", "chr_start", "chr_end"]]
    regdom = regdom[regdom["chr"].isin(list(test["chr"]))]
    for i in range(test.shape[0]):
        chrom = test.iat[i, 0]
        start = test.iat[i, 1]
        end = test.iat[i, 2]
        regdom_np = regdom["chr"].to_numpy()
        reg_start = regdom["chr_start"].to_numpy()
        reg_end = regdom["chr_end"].to_numpy()
        Chr_reduce = np.where(regdom_np == chrom)
        reg_start = np.take(reg_start, Chr_reduce, axis=0)[0]
        reg_end = np.take(reg_end, Chr_reduce, axis=0)[0]

        if any((reg_start <= start) & (reg_end >= end)):
            nb += 1
    return nb


def betacf(a, b, x):
    """Used by betai: Evaluates continued fraction for incomplete beta function"""
    maxit = 10000
    eps = 3.0e-7
    fpmin = 1.0e-30
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1
    d = 1 - qab * x / qap
    if fabs(d) < fpmin:
        d = fpmin
    d = 1 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if fabs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if fabs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if fabs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if fabs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        dell = d * c
        h *= dell
        if fabs(dell - 1.0) < eps:
            break
    if m > maxit:
        print("a or b too big, or MAXIT too small in betacf")
        return False
    return h


def betai(a, b, x):
    """Returns the incomplete beta function Ix(a, b)."""
    if x < 0 or x > 1:
        # print("bad x in routine betai")
        return False
    if x == 0 or x == 1:
        bt = 0.0
    else:
        bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1.0 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * betacf(a, b, x) / a
    return 1 - bt * betacf(b, a, 1 - x) / b


def get_binom_pval(n: int, k: int, p: float) -> float:
    """
    Calculate the binomial probability
    of obtaining k in a set of size n and whose probability is p

    Parameters
    ----------
    n : int
        Number of genomic region in the test set
    k : int
        Number of test genomic regions in the regulatory domain of a gene with annotation
    p : float
        Percentage of genome annotated

    Returns
    -------
    float
        binomial probability

    Examples
    --------
    >>> get_binom_pval(100,2,0.2)
    ...    0.9999999947037065

    """
    if k == 0:
        return 1
    else:
        return betai(k, n - k + 1, p)


def hypergeom_pmf(N: int, K: int, n: int, k: int) -> float:
    """
    Calculate the probability mass function for hypergeometric distribution

    Parameters
    ----------
    N : int
        Total number of gene in the genome
    K : int
        Number of genes in the genome with annotation
    n : int
        Number of gene in the test set
    k : int
        Number of genes in the test gene set with annotation

    Returns
    -------
    float
        proability mass function

    Examples
    --------
    >>> hypergeom_pmf(100,10,30,1)
    ...    0.11270773995748315

    """
    Achoosex = comb(K, k, exact=True)
    NAchoosenx = comb(N - K, n - k, exact=True)
    Nchoosen = comb(N, n, exact=True)
    return ((Achoosex) * NAchoosenx) / Nchoosen


def hypergeom_cdf(N: int, K: int, n: int, k: int) -> float:
    """
    Calculate the cumulative density funtion for hypergeometric distribution

    Parameters
    ----------
    N : int
        Total number of gene in the genome
    K : int
        Number of genes in the genome with annotation
    n : int
        Number of gene in the test set
    k : int
        Number of genes in the test gene set with annotation

    Returns
    -------
    float
        Cumulative density function

    Examples
    --------
    >>> hypergeom_cdf(100,10,30,1)
    ...    0.9770827595419788

    """
    return np.sum([hypergeom_pmf(N, K, n, x) for x in range(k, min(K, n) + 1)])
