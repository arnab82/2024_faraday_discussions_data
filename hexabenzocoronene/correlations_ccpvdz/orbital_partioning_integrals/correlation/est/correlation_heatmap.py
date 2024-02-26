#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Setup input arguments
parser = argparse.ArgumentParser(description='fill this out',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('infile', help='Input file to process')
args = vars(parser.parse_args())

relative = 0
edges_ref = []

fileName = args['infile']
edges = np.loadtxt(open(fileName,"rb"), delimiter=",", skiprows=0)

if relative == 1:
    edges_ref = edges
else:
    if edges_ref:
        edges = edges - edges_ref

plt.imshow(edges, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.xlabel('Column')
plt.ylabel('Row')

# plt.title('N2 Heatmap')
plt.savefig('./heatmap_correlations.png')

plt.show()  # Display the plot

# Save the plot



    