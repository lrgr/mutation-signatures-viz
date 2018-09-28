#!/usr/bin/env python

################################################################################
# CONSTANTS
################################################################################
# Load required module
import matplotlib
matplotlib.use('Agg')
import sys, os, numpy as np, seaborn as sns, matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
sns.set_style('whitegrid')

# Constants
COSMIC = 'COSMIC'
BROAD = 'Broad'
SUBS = [ 'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
N_SUBS = len(SUBS)

# Styling
SUB_COLOR = 'substitution_colors'
STYLES = {
    BROAD: {
        SUB_COLOR: {
            SUBS[0]: (122./255, 251./255., 253/255.),
            SUBS[1]: (231./255, 51./255, 35./255),
            SUBS[2]: (253./255, 252./255, 54./255),
            SUBS[3]: (147./255, 50./255, 231./255),
            SUBS[4]: (121./255, 250./255, 76./255),
            SUBS[5]: (7./255, 21./255, 244./255)
        }
    },
    COSMIC: {
        SUB_COLOR: {
            SUBS[0]: (94./255, 189./255., 235/255.),
            SUBS[1]: (5./255, 7./255, 8./255),
            SUBS[2]: (210./255, 60./255, 50./255),
            SUBS[3]: (203./255, 202./255, 203./255),
            SUBS[4]: (171./255, 205./255, 114./255),
            SUBS[5]: (230./255, 201./255, 198./255)
        }
    },
}

GREY = (203./255, 203./255, 203./255)

################################################################################
# PLOTS
################################################################################
def sbs_signature_plot(data, fig=None, sharex=False, sharey='row',
                       xlabel='Trinucleotide sequence motifs',
                       ylabel='Probability', row_labels=True,
                       title=None,
                       palette=COSMIC, fontsize=8, **kwargs):
    # Set up parameters
    if len(data.values.shape) > 1:
        K, L = data.values.shape
    else:
        K = 1
        L = len(data.values)

    if row_labels:
        row_mult = 2
        height_ratios = [ h for _ in range(K) for h in [1, 4] ]
    else:
        row_mult = 1
        height_ratios = [ h for _ in range(K) for h in [1] ]

    # Get a list of categories
    categories = data.columns
    cat_index = dict(zip(categories, range(L)))
    row_index = list(data.index)

    # Plot the number/probability of mutations per each substitution type
    fig = plt.figure(figsize=(15, K*3.1))
    gs = gridspec.GridSpec(row_mult*K, N_SUBS, height_ratios=height_ratios)
    for i in range(K):
        #  Add row labels (if necessary)
        if row_labels:
            cax = plt.subplot(gs[2*i, :])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            cax.set_facecolor(GREY)
            at = AnchoredText(row_index[i], loc=10, frameon=False, prop=dict(backgroundcolor=GREY, size=12, color='black'))
            cax.add_artist(at)

        # Add plot title (if necessary)
        if i == 0:
            if not (title is None) and not row_labels:
                plt.sup_title(title)
            elif not (title is None):
                at = AnchoredText(title, loc=6, frameon=False, prop=dict(backgroundcolor=GREY, size=12, weight='bold', color='black'))
                cax.add_artist(at)

        # Plot distributiton per substitution type
        for j, sub in enumerate(SUBS):
            # Get the list of active categories
            sub_cats = [ c for c in categories if sub in c ]
            sub_cat_idx = [ cat_index[c] for c in sub_cats ]
            xticklabels = [ c.replace('[%s]' % sub, '-') for c in sub_cats ]

            # Compute the number of mutations across this substitution category
            x = np.arange(len(sub_cats))
            y = data.values[i, sub_cat_idx]

            # Plot
            if (sharey == 'row' and j == 0) or (sharey and j == 0 and i == 0):
                first_ax = ax = plt.subplot(gs[2*i+1, j])
            else:
                ax = plt.subplot(gs[2*i+1, j], sharey=first_ax)

            ax.bar(x, y, align='center', color=STYLES[palette][SUB_COLOR][sub])
            ax.set_title(sub)

            # Only show at most 16 xticklabels
            if len(sub_cats) > 16:
                offset = int(round(len(sub_cats)/16))
                ax.grid(False)
                ax.set_xticks(x[::offset])
                ax.set_xticklabels(xticklabels[::offset], rotation='vertical', fontsize=fontsize)
            else:
                ax.set_xticks(x)
                ax.set_xticklabels(xticklabels, rotation='vertical', fontsize=fontsize)

            #  Add axis labels and ticks (sharing as necessary)
            if i == K-1:
                ax.set_xlabel(xlabel)

            if j == 0:
                ax.set_ylabel(ylabel)
            elif sharey == 'row' or sharey:
                plt.setp(ax.get_yticklabels(), visible=False)
