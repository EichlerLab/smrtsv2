# Code based on:
# http://nbviewer.ipython.org/github/rasbt/pattern_classification/blob/master/clustering/hierarchical/clust_complete_linkage.ipynb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
import argparse
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
import sys

sys.setrecursionlimit(10000)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genotypes")
    parser.add_argument("output")
    args = parser.parse_args()

    df = pd.read_table(args.genotypes).T

    hierarchy.set_link_color_palette(['black'])
    fig = plt.figure(figsize=(10, 6))

    # Row clusters
    axd1 = fig.add_axes([0.0, 0.1, 0.07, 0.6])

    print "Row distance"
    row_dists = pdist(df, metric='euclidean')
    print "Row linkage"
    row_clusters = linkage(df, method='average')
    print "Row dendrogram"
    row_dendr = dendrogram(
        row_clusters,
        distance_sort='ascending',
        orientation='right',
        color_threshold=np.inf
    )

    axd1.set_xticks([])
    axd1.set_yticks([])

    # Column clusters
    axd2 = fig.add_axes([0.26, 0.74, 0.6, 0.20])

    print "Column clusters"
    col_dists = pdist(df.T, metric='euclidean')
    print "Column linkage"
    col_clusters = linkage(col_dists, method='single')
    print "Column dendrogram"
    col_dendr = dendrogram(
        col_clusters,
        orientation='top',
        no_plot=True
    )

    axd2.set_xticks([])
    axd2.set_yticks([])

    print "Reorder rows and columns"
    df_rowclust = df.ix[row_dendr['leaves'][::-1]]
    df_rowclust.columns = [df_rowclust.columns[col_dendr['leaves']]]

    # Remove axes spines from dendrogram
    for i,j in zip(axd1.spines.values(), axd2.spines.values()):
        i.set_visible(False)
        j.set_visible(False)

    axm = fig.add_axes([0.26,0.1,0.6,0.6]) # x-pos, y-pos, width, height
    for region in ("top", "right", "bottom", "left"):
        axm.spines[region].set_visible(False)

    cax = axm.matshow(df_rowclust, interpolation="nearest", cmap=pylab.cm.Blues, vmin=0, vmax=1, aspect="auto")

    fig.canvas.draw()
    labels = [item.get_text() for item in axm.get_xticklabels(minor=True)]
    print labels
    print df_rowclust.index
    axm.set_xticklabels(labels, size=8, family="serif")

    y_axis_labels = ["%s (%s)" % i for i in zip(list(df_rowclust.index), df_rowclust.sum(axis=1))]
    axm.yaxis.set_ticks(np.arange(0, len(df_rowclust.index), 1))
    axm.set_yticklabels(y_axis_labels, size=8, family="serif")

    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        left='off',
        right='off',
        labelbottom='off'
    )

    #plt.tight_layout()
    plt.savefig(args.output)
