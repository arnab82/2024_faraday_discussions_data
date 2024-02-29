import numpy as np
import matplotlib.pyplot as plt

matrix = np.array([[-2.071532588, -2.065048656, -2.06626555],
                   [-2.463916407, -2.489055602, -2.490800125],
                   [-2.473229996, -2.493793923, -2.495326741]])

num_rows, num_cols = matrix.shape

x = np.arange(2, num_rows+2,1)
width = 0.2

fig, ax = plt.subplots()
ax.set_ylim(0, -3.5)


for i in range(num_cols):
    ax.bar(x + i*width, matrix[:, i], width, label=f'${{\delta_{{elec}}}}={i+1}$')


ax.set_ylabel('Correlation Energy (eV)')
ax.set_xlabel("nbody")
# ax.set_title("Contribution of nbody terms with the ${\delta e^{-1}}$ to the correlation energy")
ax.set_xticks(x)
ax.legend(fontsize=12, loc='upper left')

# Add bar labels
for i in range(num_rows):
    for j in range(num_cols):
        ax.text(x[i] + (j+0.5)*width, matrix[i, j], f'{matrix[i, j]:.4f}', ha='right', va='bottom',rotation=90)
fig.set_dpi(300)
fig.set_size_inches(5,5)
fig.savefig('nbody_correlations.jpg', bbox_inches='tight', dpi=300)
plt.show()






