#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt


relative = 0
edges_ref = []
# edges = np.array([[-0.0811842, -0.08107293, -0.08109712],
#                    [-0.10069981, -0.10076323, -0.1008057],
#                    [-0.1013064, -0.10128546, -0.10131803]])
edges=np.array([[-47.77015933,	-47.62063792,	-47.64869985],
[-56.81864722,	-57.39836456,	-57.43859378],
[-57.03342137,	-57.50763166,	-57.5429789]])
edges=np.array([[-2.071532588,	-2.065048656,	-2.06626555],
[-2.463916407,	-2.489055602,	-2.490800125],
[-2.473229996,	-2.493793923,	-2.495326741]])
num_rows, num_cols = edges.shape
if relative == 1:
    edges_ref = edges
else:
    if edges_ref:
        edges = edges - edges_ref
    plt.imshow(edges, cmap='RdGy', interpolation='nearest', vmin=-2.5, vmax=-2.05, aspect='auto')
    plt.colorbar()      
    plt.xlabel('Nbody')
    plt.ylabel("${\delta}$-electron")
    plt.xticks([0, 1, 2], ['2', '3', '4'])
    plt.yticks([0, 1, 2], ['1', '2', '3'])
    plt.savefig('./heatmap_nbody.png')
    plt.show()  

 




    