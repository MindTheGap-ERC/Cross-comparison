#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from comparePlot import floatMarl

folders = [
    'data/larger_system/Fortran/Tstar_100_L_625',
    'data/larger_system/Fortran/Tstar_100_L_750',
    'data/larger_system/Fortran/Tstar_100_L_875',
    'data/larger_system/Fortran/Tstar_100_L_1000',
    'data/larger_system/Fortran/Tstar_100_L_1125',
    'data/larger_system/Fortran/Tstar_5000_L_1625',
]

labels = ["Aragonite", "Calcite", "Porosity", r"$Ca^{2+}$", r"$CO_3^{2-}$"]
line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
colors = ['#E69F00', 'red', '#000000', '#009E73', '#0072B2']

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

depth_ind = 4

Ts = 131.9 / 0.1**2

# Plot for each folder
for idx, folder in enumerate(folders):
    filename = f'{folder}/amarlx'
    ax = axes[idx]

    df = pd.read_csv(filename, sep=r'\s+')
    floatMarl(df)
    t_plot = df.t * Ts / 1000

    ax.plot(t_plot, df[df.columns[depth_ind]], label=labels[0],
            linestyle=line_styles[0], color=colors[0])
    ax.plot(t_plot, df[df.columns[depth_ind+4]], label=labels[1],
            linestyle=line_styles[1], color=colors[1])
    ax.plot(t_plot, df[df.columns[depth_ind+8]], label=labels[2],
            linestyle=line_styles[2], color=colors[2])
    ax.plot(t_plot, df[df.columns[depth_ind+12]], label=labels[3],
            linestyle=line_styles[3], color=colors[3])
    ax.plot(t_plot, df[df.columns[depth_ind+16]], label=labels[4],
            linestyle=line_styles[4], color=colors[4])

    # Extract parameters from folder name
    folder_name = folder.split('/')[-1]
    parts = folder_name.split('_')
    tstar = parts[1]
    L = parts[3]

    ax.set_title(f'$T^*$ = {tstar}, L = {L}', fontsize=10)
    ax.set_xlabel('Time [ky]')
    ax.set_ylabel('Concentration/Porosity')

    # panel label in upper left corner
    panel_label = chr(97 + idx)  
    ax.text(0.02, 0.98, f'({panel_label})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')

    if idx == 5:
        ax.legend(loc='lower right', fontsize=8)

plt.tight_layout()
plt.savefig('figs/larger_system/Fig_A1_larger_system.png',
            format='png', dpi=300, bbox_inches='tight')

plt.close()
