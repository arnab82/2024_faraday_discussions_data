import os
import sys
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors

import scipy

file = ""
if len(sys.argv) > 1:
    file = os.path.splitext(sys.argv[1])[0]
else:
    print("No input file provided.")

print(" Working with file: ", file)


np.set_printoptions(suppress=True, precision=6, linewidth=1500)

data_path = file+'.csv'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader)
    data = list(reader)

conversion = 219474.63 # cm-1
conversion = 1 # hartrees
conversion = 1000 # mH 

energy_var = {}
energy_pt2 = {}
num_thresh = 6
       
assert(len(data[0]) % 2 == 0)

n_extrap_points = int((len(data[0]))//2)
print(" Number of extrapolation points: ", n_extrap_points)

for i in range(len(data)):
    if  i > 0:
        energy_var[i] = np.array([float(a)*conversion for a in data[i][0:n_extrap_points]])
        energy_pt2[i] = np.array([float(a)*conversion for a in data[i][n_extrap_points:2*n_extrap_points]])

print(energy_var)
print(energy_pt2)


blue = '#3E6D9C'
orange = '#FD841F'
red_orange = '#E14D2A'
dark_blue = '#001253'

cb = ['#000000']
cb.extend([i for i in plt.rcParams['axes.prop_cycle'].by_key()['color']])

fig, ax = plt.subplots()

#this sets the number of decimal points on axis energies
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
extrap = []


#set your x and y axis limits

xmax = 0 
xmin = -20

ymax = 30
ymin = -0.0



# Find shift for plots to have the extrapolation at the origin
shift = 0.0
for key in energy_var:
    x = energy_pt2[key] - energy_var[key]
    y = energy_var[key]
    m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)

    shift = min(b, shift)

for i in energy_pt2:
    energy_pt2[i] -= shift
    energy_var[i] -= shift 



for key in energy_var:
    x = energy_pt2[key] - energy_var[key]
    z = energy_pt2[key]
    y = energy_var[key]
    

    m, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    m2, b2, r_value, p_value, std_err = scipy.stats.linregress(x, z)
    
    extrap.append(b)

    plt.rcParams.update({'font.size': 10})

    #ymin = min(ymin, b)
    #ymax = max(ymax, m*xmin+b)

    ax.plot(x, y, marker='o',linestyle='' ,markersize=8, color = cb[key])
    ax.plot(x, z, marker='x',linestyle='' , markersize=8, color = cb[key])

    x2 = np.array([-1,0])*conversion
    print(x2)
    line = m*x2+b
    print(line)
    ax.plot(x2, line, alpha=1.0, color = cb[key], linestyle='-', linewidth=1.5,label='TPSCI' if key == 1 else 'HCI')
    line = m2*x2+b                                     
    ax.plot(x2, line, alpha=0.5, color = cb[key], linestyle='--', linewidth=1.5,label='TPSCI+PT2' if key == 1 else 'HCI+PT2')
    print("Extrapolated Result: %14.8f"% ((b+shift)/conversion))
    print("R^2                : %14.8f"% r_value)
    print("Var root",key,y)
    print("PT  root",key,z)
    #print("DIFFF",x)
    #print("color",cb[key])


ymin = ymin 
print("x: ", (xmin, xmax))
print("y: ", (ymin, ymax))
ax.set_ylim(ymin, ymax)
ax.set_xlim(xmin, xmax)

ax.set_xlabel('$\Delta$E$_{PT2}$ (mH ) ')
ax.set_ylabel('Energy (mH) ')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
# ax.set_yticklabels([])
#ax.legend()
ax.set_title('Extrapolated TPSCI Energy, %12.8f (Hartree) '%(shift/conversion), fontsize=10)
#ax.set_title("Tetracence Tetramere \n (40o, 40e) \n ε = {0.0005, 0.0007, 0.001, 0.005, 0.01} \n 31 roots")
# blue_patch = mpatches.Patch(color=blue, label='TPSCI')
# blue_dotted = mpatches.Patch(color=blue,linestyle='--', label='TPSCI+PT2')
# red_patch = mpatches.Patch(color=red_orange, label='HCI')
# red_dotted = mpatches.Patch(color=red_orange, label='HCI+PT2', linestyle='--')
# ax.legend(handles=[blue_patch,blue_dotted,red_patch,red_dotted], loc='upper right')
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()

ax.legend()
fig = plt.gcf()
fig.set_size_inches(4.5,4.5)
fig.savefig(file+'_TPSCIextrap_ccpvdz.pdf', dpi=300, bbox_inches='tight')
fig.savefig(file+'_TPSCIextrap_ccpvdz.png', dpi=300, bbox_inches='tight')

plt.show()

print(extrap)
