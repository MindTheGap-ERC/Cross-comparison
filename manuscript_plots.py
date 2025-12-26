import matplotlib.pyplot as plt
import pandas as pd
from comparePlot import floatMarl

# Create figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Common settings
depth_ind = 4
Ts = 131.9/0.1**2  # time scaling constant
colors = ['#E69F00', 'red', '#000000', '#009E73', '#0072B2']
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
labels = ["Aragonite", "Calcite", "Porosity", r"$Ca^{2+}$", r"$CO_3^{2-}$"]

# Plot first dataset (left panel - Fig 4a)
df1 = pd.read_csv("data/replication/Fortran/Fig._4a/amarlx", sep=r'\s+')
floatMarl(df1)
t_plot1 = df1.t * Ts / 1000

for i, (label, ls, color) in enumerate(zip(labels, linestyles, colors)):
    ax1.plot(t_plot1, df1[df1.columns[depth_ind + i*4]],
            label=label, linestyle=ls, color=color)

ax1.set_xlabel('Time [ky]')
ax1.set_ylabel('Concentration/Porosity')
ax1.text(0.02, 0.98, '(a)', transform=ax1.transAxes,
        fontsize=12, fontweight='bold', va='top')

# Plot second dataset (right panel - Fig 4b)
df2 = pd.read_csv("data/replication/Fortran/Fig._4b/amarlx", sep=r'\s+')
floatMarl(df2)
t_plot2 = df2.t * Ts / 1000

for i, (label, ls, color) in enumerate(zip(labels, linestyles, colors)):
    ax2.plot(t_plot2, df2[df2.columns[depth_ind + i*4]],
            label=label, linestyle=ls, color=color)

ax2.set_xlabel('Time [ky]')
ax2.set_ylabel('Concentration/Porosity')
ax2.legend(loc='lower right')
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes,
        fontsize=12, fontweight='bold', va='top')

# Adjust layout and save
plt.tight_layout()
plt.savefig('Fortran_temporal_panel.svg', format='svg', bbox_inches='tight')
plt.close()
