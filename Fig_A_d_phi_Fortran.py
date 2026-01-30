#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from comparePlot import floatMarl

# Folders ordered by d_phi value
folders = [
    'data/d_phi/Fortran/0.001',
    'data/d_phi/Fortran/0.005',
    'data/d_phi/Fortran/0.01',
    'data/d_phi/Fortran/0.02',
    'data/d_phi/Fortran/0.03',
    'data/d_phi/Fortran/0.04',
    'data/d_phi/Fortran/0.05',
    'data/d_phi/Fortran/0.1',
]

# Common settings from FORTRAN RUNS section
depth_ind = 4
Ts = 131.9 / 0.1**2  # time scaling constant
colors = ['#E69F00', 'red', '#000000', '#009E73', '#0072B2']
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
labels = ["Aragonite", "Calcite", "Porosity", r"$Ca^{2+}$", r"$CO_3^{2-}$"]

# Create figure with 2x4 subplots
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()

# Collect all y-values to determine common limits
all_y_values = []

for idx, folder in enumerate(folders):
    ax = axes[idx]

    # Read the amarlx file
    filename = f"{folder}/amarlx"
    df = pd.read_csv(filename, sep=r'\s+')
    floatMarl(df)

    t_plot = df.t * Ts / 1000

    # Plot all variables
    for i, (label, ls, color) in enumerate(zip(labels, linestyles, colors)):
        y_data = df[df.columns[depth_ind + i*4]]
        ax.plot(t_plot, y_data, label=label, linestyle=ls, color=color)
        all_y_values.extend(y_data)

    # Extract d_phi value from folder name
    d_phi = folder.split('/')[-1]

    ax.set_title(rf'$D_\phi$ = {d_phi}', fontsize=10)
    ax.set_xlabel('Time [ky]')
    ax.set_ylabel('Concentration/Porosity')

    # Add panel label
    panel_label = chr(97 + idx)
    ax.text(0.02, 0.98, f'({panel_label})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')

    # Add legend to last panel
    if idx == 7:
        ax.legend(loc='lower right', fontsize=8)

# Set common y-axis limits
y_min = min(all_y_values)
y_max = max(all_y_values)
y_margin = (y_max - y_min) * 0.05
for i in range(8):
    axes[i].set_ylim(y_min - y_margin, y_max + y_margin)

plt.tight_layout()
plt.savefig('figs/Fig_A_d_phi_Fortran.png',
            format='png', dpi=300, bbox_inches='tight')

plt.close()
