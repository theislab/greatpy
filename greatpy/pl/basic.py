import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors, rcParams
from numpy import log
import matplotlib
from matplotlib import rc 

import greatpy as gp

plt.rcParams.update({"font.size": 14, "font.weight": "normal"})

def scatterplot(
    great_df: pd.DataFrame, colname_x: str, colname_y: str,
    title: str = "", minus_log10=True, ax=None
) -> None:
    """
    Create a scatterplot from a
    pandas dataframe between two columns.
    A logarithmic scale can be used.

    Parameters
    ----------
    great_df : pd.DataFrame
        Output of the greatpy.tl.GREAT.enrichment function
    colname_x : str
        Name of the column to be used as x axis
    colname_y : str
        Name of the column to be used as y axis
    title : str
        Title of the plot
    minus_log10 : bool
        If True, the logarithmic scale is used
    ax :
        Define the position of the plot in a figure

    Returns
    -------
    None

    """
    great_df = great_df.dropna()
    great_df = great_df.astype({colname_x: "float64", colname_y: "float64"})
    if minus_log10:
        great_df[f"-log({colname_x})"] = -log(great_df[colname_x])
        great_df[f"-log({colname_y})"] = -log(great_df[colname_y])
        sns.scatterplot(data=great_df, x=f"-log({colname_x})", y=f"-log({colname_y})", ax=ax).set_title(title)
    else:
        sns.scatterplot(data=great_df, x=colname_x, y=colname_y, ax=ax).set_title(title)

def graph_nb_asso_per_peaks(
    test: str or pd.DataFrame, regdom: str or pd.DataFrame,
    ax=None, color=None) -> None:
    """
    Creates a barplot representing the
    percentage of peaks for all possible association numbers

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested
    regdom : str or pd.DataFrame
        Regulatory domain of all genes in the genome
    ax :
        Define the position of the plot in a figure

    Returns
    -------
    None

    """
    nb_asso_per_peaks = gp.tl.get_nb_asso_per_region(test, regdom)

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
    g.set_xlabel("Number of associated genes per region", fontsize=13)
    g.set_ylabel("Genomic region (%)", fontsize=13)

    for i in range(nb.shape[0]):
        x = nb.iloc[i]["number"]
        y = nb.iloc[i]["percentage"]
        g.text(x=x - 0.06, y=y + 1, s=nb.number_genes[i])

def graph_dist_tss(
    test: str or pd.DataFrame, regdom: str or pd.DataFrame, 
    ax=None, color="#325fa8") -> None:
    """
    Creation of a barplot of the distance
    between the peaks and the TSS of the associated gene(s).

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested
    regdom : str or pd.DataFrame
        Regulatory domain of all genes in the genome
    ax :
        Define the position of the plot in a figure

    Returns
    -------
    None

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

    dist = gp.tl.get_dist_to_tss(test, regdom)
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
    g.set_xlabel("Distance to TSS (kb)", fontsize=13)
    g.set_ylabel("Genomic region (%)", fontsize=13)
    g.set_title("Binned by absolute distance to TSS", fontsize=20)

def graph_absolute_dist_tss(
    test: str or pd.DataFrame, regdom: str or pd.DataFrame, 
    ax=None, color="#325fa8") -> None:
    """
    Creation of a barplot of the absolute
    distance between the peaks and the TSS of the associated gene(s).

    Parameters
    ----------
    test : str or pd.DataFrame
        Genomic set of peaks to be tested
    regdom : str or pd.DataFrame
        Regulatory domain of all genes in the genome
    ax :
        Define the position of the plot in a figure

    Returns
    -------
    None

    """
    res = {"0:5": [0], "5:50": [0], "50:500": [0], ">500": [0]}
    nb = 0
    dist = gp.tl.get_dist_to_tss(test, regdom)

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
    g.set_xlabel("Absolute distance to TSS (kb)", fontsize=13)
    g.set_ylabel("Genomic region (%)", fontsize=13)
    g.set_title("Binned by absolute distance to TSS", fontsize=20)

def scale_data_5_75(data):
    mind = np.min(data)
    maxd = np.max(data)

    if maxd == mind:
        maxd = maxd + 1
        mind = mind - 1

    drange = maxd - mind
    return (((data - mind) / drange * 0.70) + 0.05) * 100

def plot_enrich(data, n_terms=20, color="cool", save=False):
    """
    Creation of a dotplot of the enrichment
    GO term in the inputs datas

    Parameters
    ----------
    data : pd.DataFrame
        Results for greatpy
    n_terms : int
        the number of term to be shown
    color : str
        The color of the cmap in the plot
    save : bool
        Is the plot should be save

    Returns
    -------
    None

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

    rcParams.update({"font.size": 14, "font.weight": "bold"})

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
    ax.set_xlabel("Gene ratio", fontsize=14, fontweight="bold")
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
    cbar.set_label("Adjusted p-value", fontsize=14, fontweight="bold")

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
        plt.savefig("dotplot_save", dpi=500)

    plt.show()

def make_bubble_heatmap(order_frame, sizeDict, na_color='gray', title='title',
                        tickscolorbar=[-2, -1, 0, 1, 2], vmin=-2.5, vmax=2.5,
                        heatmap_grid=[2, 4, 0, 2, 2, 1],
                        circle_legend_grid=[2, 4, 0, 2, 2, 1],
                        colorbar_grid=[2, 5, 0, 3, 2, 1],
                        palette_id='RdBu_r', cbar_label='cbar_label', ncols=8,
                        marker=None, **kwargs):
    """
    Generate a dotplot with multiple categories

    Parameters
    ----------
    order_frame
        DataFrame with gene ontologies 
    sizeDict
        Dictionary with gene ontologies as keys and sizes as values
    na_color
        Color for NA values
    title
        Title for the plot
    tickscolorbar
        Tick marks for the colorbar
    vmin
        Minimum value for the colorbar
    vmax
        Maximum value for the colorbar
    heatmap_grid
        Grid for the heatmap
    circle_legend_grid
        Grid for the circle legend
    colorbar_grid
        Grid for the colorbar
    palette_id
        Palette id for the colorbar
    cbar_label
        Label for the colorbar
    ncols
        Number of columns for the colorbar
    marker
        Marker for the dots
    kwargs
        Additional keyword arguments for the plot

    Returns
    -------
    None 
        Bubble heatmap plot for gene ontology
    """
    from matplotlib import gridspec
    sns.set_style('white')

    rc("font",weight="normal")
    # rc("alpha",alpha=0.8)

    quantAmplifier = kwargs.get('quantAmplifier', 1) #factor for size of bubbles

    #function you received on the color gradient
    pallete = None
    if palette_id == 'RdBu_r':
        pallete = plt.cm.RdBu_r
    elif palette_id == 'Blues':
        pallete = plt.cm.Blues
    elif palette_id == 'Reds':
        pallete = plt.cm.Reds
    elif palette_id == 'Greens':
        pallete = plt.cm.Greens
    elif palette_id == 'Purples':
        pallete = plt.cm.Purples
    elif palette_id == 'YlGnBu':
        pallete = plt.cm.YlGnBu
    elif palette_id == 'PuOr':
        pallete = plt.cm.PuOr
    else:
        palette = palette_id
    scalarmap,colorList = get_specific_color_gradient(pallete,
        np.array(order_frame),
        vmin=vmin,
        vmax=vmax)

    fig = kwargs.get('fig', None)
    if fig is None:
        plt.clf()
        fig = plt.figure(figsize=(kwargs.get('w', 11), kwargs.get('h', 5)))


    nrows, ncols, rowi, coli, rowspan, colspan = heatmap_grid
    gs = gridspec.GridSpec(nrows, ncols)
    ax = plt.subplot(gs[rowi: rowi + rowspan, coli: coli + colspan])

    ylabelList = []
    i = 0

    # to keep the order of dataframes in the same order as sizes, revert the dataframes
    order_frame = order_frame.reindex(index=order_frame.index[::-1])
    sizeDict = sizeDict.reindex(index=sizeDict.index[::-1])

    min_circle_size = kwargs.get('min_circle_size')
    max_circle_size = kwargs.get('max_circle_size')

    for ri, r in order_frame.iterrows():
        patList = list(r.values)

        sizeList = []
        if min_circle_size is not None and max_circle_size is not None:
            for si in list(sizeDict.loc[ri]):
                if si < min_circle_size:
                    si = min_circle_size
                scaled = (si - min_circle_size) / (max_circle_size - min_circle_size)
                sizeList.append(scaled * quantAmplifier)
        else:
            sizeList = [(si ** kwargs.get('power', 1.0)) * (quantAmplifier) for si in list(sizeDict.loc[ri])]

        colorList = scalarmap.to_rgba(patList)
        colorList = [ci if not np.isnan(vi) else na_color
                        for ci, vi in zip(colorList, patList)]
        edgecolorList = list()
        x = list(range(len(patList)))
        y = [i] * len(patList)

        # define hataches, linewidths and alphas
        hatches = [None if np.isnan(pi) else None for pi in patList]
        linewidths = [0.1 if np.isnan(pi) else kwargs.get('sig_line_width', 2.0) for pi in patList]\
            if kwargs.get('line_widths', None) is None else list(kwargs.get('line_widths').loc[ri])
        alphas = [.5 if np.isnan(alpha) else 1.0 for alpha in patList]

        for xi, yi, si, ci, hatch_i,lw, alpha in zip(x, y, sizeList, colorList,hatches, linewidths, alphas):
            ax.scatter(xi, yi, s=abs(si) * quantAmplifier,
                        marker=marker if marker is not None else ('v' if si < 0 else '^'),
                        hatch=hatch_i, alpha=alpha,
                        edgecolor=kwargs.get('edgecolor', 'black'), color=ci, linewidth=lw)
        ylabelList.append(ri)
        i += 1
    if kwargs.get('grid', True):
        plt.grid(True, linewidth=kwargs.get('grid_linewidth', 1.0))
    plt.title(kwargs.get('heatmap_title', 'title'), fontname='Arial',fontweight='normal')
    plt.xlabel(kwargs.get('xlab', 'xlab'), fontsize=14,fontname='Arial',fontweight='normal')
    plt.ylabel(kwargs.get('ylab', 'ylab'), fontsize=14,fontname='Arial',fontweight='normal')


    plt.yticks(list(range(len(ylabelList))))
    ax.set_yticklabels(ylabelList, fontsize=kwargs.get('yticks_fontsize', 11), color="black",
                        ha=kwargs.get('ha_ylabs', 'right'),fontname='Arial',fontweight='normal',alpha = 1)

    # # set different color for x_tick labels 
    # nb = int((len(ylabelList)+1)/(len(order_frame.columns)+1))
    # color = (["r"]*5+["black"]*5)*((len(order_frame.columns)+1) // 2)
    # if len(color) != len(order_frame.columns): 
    #     color += ["r"]*(len(order_frame.columns)-len(color))
    
    # for ticklabel,tickcolor in zip(plt.gca().get_yticklabels(),color):
    #     ticklabel.set_color(tickcolor)

    remove_top_n_right_ticks(ax)

    plt.xticks(list(range(len(order_frame.columns))), order_frame.columns,
                fontsize=kwargs.get('xticks_fontsize', 11), rotation=kwargs.get('rotation_xlabs', 90), color="black",
                ha=kwargs.get('ha_xlabs', 'center'),fontname='Arial',fontweight='normal')

    lh, lt = get_legendHandle_for_second_sanity_check_plot(quantAmplifier=quantAmplifier * quantAmplifier,
                                                            marker='^' if marker is None else marker,
                                                            fmt=kwargs.get('fmt_legend', '%.1f'),
                                                            lw=0.5,
                                                            min_circle_size=min_circle_size,
                                                            max_circle_size=max_circle_size,
                                                            values=[tick ** kwargs.get('power', 1.0) for tick in kwargs.get('circle_legend_ticks')],
                                                            labels=kwargs.get('circle_legend_ticks'))
    if kwargs.get('show_circle_legend', True):
        l = plt.legend(lh, lt, bbox_to_anchor=kwargs.get('circle_legend_bbox', (1.8, 1)), scatterpoints=1,
                        title=kwargs.get('circles_legend_title', 'circles_legend_title'), ncol=1,
                        frameon=False)
        # Add the legend manually to the current Axes.
        ax = plt.gca().add_artist(l)

    # this is to add the circle (significant or not
    if kwargs.get('show_sig_legend', False):
        lh, lt = get_legendHandle_for_second_sanity_check_plot(quantAmplifier=quantAmplifier * quantAmplifier,
                                                                marker='^' if marker is None else marker,
                                                                fmt=kwargs.get('sig_fmt_legend', '%s'),
                                                                lw=[kwargs.get('sig_line_width', 2.0), 0.0],
                                                                labels=kwargs.get('sig_legend_ticks',
                                                                                    ['Yes', 'No']),
                                                                values=kwargs.get('sig_legend_ticks',
                                                                                    ['Yes', 'No']),
                                                                edgecolor=kwargs.get('edgecolor', 'black'),
                                                                min_size_default=max(kwargs.get('circle_legend_ticks')) ** kwargs.get('power', 1.0))

        plt.legend(lh, lt, bbox_to_anchor=kwargs.get('sig_legend_bbox', (2.4, 1)),
                    title=kwargs.get('sig_legend_title', 'significant'), ncol=1, scatterpoints=1,
                    frameon=False)
    
    # nrows, ncols, rowi, coli, rowspan, colspan = colorbar_grid
    # plt.show()
    # ax1 = plt.subplot2grid([nrows, ncols], [rowi, coli], rowspan=rowspan, colspan=colspan)
    
    # plt.axis('off')
    # ax1.set_xticklabels([])
    # ax1.set_yticklabels([])

    if kwargs.get('show_colorbar', True):
        cbar = fig.colorbar(scalarmap, orientation="horizontal", format=kwargs.get('cbar_fmt_ticks', "%.1f"),
                            ticks = tickscolorbar,pad=0.05)
        cbar.ax.tick_params(labelsize=kwargs.get('colorbar_ticks_labelsize', 12))
        cbar.set_label(cbar_label, fontsize=12,fontweight='normal')
    despine_all()

def get_specific_color_gradient(colormap,inputList,**kwargs):
    vmin = kwargs.get('vmin','blaq')
    vmax = kwargs.get('vmax','blaq')
    cm = plt.get_cmap(colormap)
    if vmin=='blaq' or vmax=='blaq':
        if type(inputList)==list:
            cNorm = matplotlib.colors.Normalize(vmin=min(inputList), vmax=max(inputList))
        else:
            cNorm = matplotlib.colors.Normalize(vmin=inputList.min(), vmax=inputList.max())
    else:
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax = vmax)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    scalarMap.set_array(inputList)
    colorList=scalarMap.to_rgba(inputList)
    return scalarMap,colorList

def despine_all():
    sns.despine(offset=10, trim=True, top=True, right=True, left=True,
                bottom=True)

def remove_top_n_right_ticks(ax):
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

def get_legendHandle_for_second_sanity_check_plot(quantAmplifier=None,
                                                  labels=None,
                                                  values=None,
                                                  marker='^',
                                                  fmt='%.1f',
                                                  lw=2.0, min_circle_size=None, max_circle_size=None,
                                                  min_size_default=1.0,
                                                  edgecolor='black',
                                                  color='grey'):
    labels = [0.1, 0.5, 1.0, 1.5] if labels is None else labels

    legendHandleList = list()
    labelSize = list()
    for i, quantMem in enumerate(values):
        v = values[i]

        if min_circle_size is not None and max_circle_size is not None:
            scaled = (quantMem - min_circle_size) / (max_circle_size - min_circle_size)
            quantMem = scaled

        next_size = (quantMem * quantAmplifier) if not isinstance(quantMem, str) else min_size_default * quantAmplifier

        legendHandleList.append(plt.scatter([], [], s=next_size,
                                            color=color, edgecolor=edgecolor,
                                            linewidths=lw if not isinstance(lw, list) else lw[i],
                                            alpha=0.9, marker=marker))
        labelSize.append(fmt % labels[i])
    return legendHandleList, labelSize