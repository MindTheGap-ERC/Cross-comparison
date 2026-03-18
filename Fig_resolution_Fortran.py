#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from comparePlot import floatMarl

folders = [
    'data/resolution/FORTRAN/Tstar_30_nnx_200',
    'data/resolution/FORTRAN/Tstar_10_nnx_300',
    'data/resolution/FORTRAN/Tstar_10_nnx_400',
]

depth_ind = 4  # base of domain
Ts = 131.9 / 0.1**2  # time scaling constant
t_max_ky = 10 * Ts / 1000  # 10 T* in ky, for truncating the Tstar_30 run

colors = ['#E69F00', 'red', '#000000', '#009E73', '#0072B2']
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
labels = ["Aragonite", "Calcite", "Porosity", r"$Ca^{2+}$", r"$CO_3^{2-}$"]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

all_y_values = []

for idx, folder in enumerate(folders):
    ax = axes[idx]

    df = pd.read_csv(f'{folder}/amarlx', sep=r'\s+')
    floatMarl(df)

    t_plot = df.t * Ts / 1000

    # Truncate Tstar_30 run to first 10 T*
    mask = t_plot <= t_max_ky
    t_plot = t_plot[mask]
    df = df[mask]

    for i, (label, ls, color) in enumerate(zip(labels, linestyles, colors)):
        y_data = df[df.columns[depth_ind + i * 4]]
        ax.plot(t_plot, y_data, label=label, linestyle=ls, color=color)
        all_y_values.extend(y_data)

    folder_name = folder.split('/')[-1]
    parts = folder_name.split('_')
    nnx = parts[3]

    ax.set_title(f'nnx = {nnx}', fontsize=10)
    ax.set_xlabel('Time [ky]')
    ax.set_ylabel('Concentration/Porosity')

    panel_label = chr(97 + idx)
    ax.text(0.02, 0.98, f'({panel_label})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')

    if idx == 2:
        ax.legend(loc='lower right', fontsize=8)

y_min = min(all_y_values)
y_max = max(all_y_values)
y_margin = (y_max - y_min) * 0.05
for ax in axes:
    ax.set_ylim(y_min - y_margin, y_max + y_margin)

plt.tight_layout()
plt.savefig('figs/Fig_resolution_Fortran.png',
            format='png', dpi=300, bbox_inches='tight')

plt.close()
