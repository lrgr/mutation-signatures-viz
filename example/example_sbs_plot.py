#!/usr/bin/env Python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.join(this_dir, '../src'))
from mutation_signatures_visualization import sbs_signature_plot

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-if', '--input_file', type=str, required=True)
parser.add_argument('-of', '--output_file', type=str, required=True)
parser.add_argument('-t', '--title', type=str, required=True)
parser.add_argument('-as', '--active_signatures', type=int, required=False,
                    nargs='*', default=[], help='1-based indices.')
args = parser.parse_args(sys.argv[1:])

# Load the signatures
sig_df = pd.read_csv(args.input_file, sep='\t', index_col=0)
K, L = sig_df.values.shape
print('Loaded %s signatures across %s mutation categories' % (K, L))

# Restrict to "active" signatures
if len(args.active_signatures) > 0:
    sig_df = sig_df.loc[[ sig_df.index[idx-1] for idx in args.active_signatures ]]
    print('Restricted to signatures: %s' % ', '.join(map(str, args.active_signatures)))

# Plot and save to file
sbs_signature_plot(sig_df, title=args.title, sharey=True)
plt.tight_layout(pad=2)
plt.savefig(args.output_file)
