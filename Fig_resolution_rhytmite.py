#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob

folders = [
    'data/resolution/rhythmite/200nnx',
    'data/resolution/rhythmite/300nnx',
    'data/resolution/rhythmite/400nnx',
]

labels = ["Aragonite", "Calcite", "Porosity", r"$Ca^{2+}$", r"$CO_3^{2-}$"]
line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
colors = ['#E69F00', 'red', '#000000', '#009E73', '#0072B2']

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

Ts = 131.9 / 0.1**2

all_y_values = []

for idx, folder in enumerate(folders):
    ascii_files = glob.glob(f'{folder}/*.ascii')
    if not ascii_files:
        print(f"Warning: No .ascii file found in {folder}")
        continue
    filename = ascii_files[0]
    ax = axes[idx]

    df = (pd.read_csv(filename, sep=r'\s+')).shift(axis=1).iloc[:,1:]
    t_plot = np.array(df.x * Ts / 1000)

    ax.plot(t_plot, np.array(df.AR),  label=labels[0], color=colors[0], linestyle=line_styles[0])
    ax.plot(t_plot, np.array(df.CA),  label=labels[1], color=colors[1], linestyle=line_styles[1])
    ax.plot(t_plot, np.array(df.phi), label=labels[2], color=colors[2], linestyle=line_styles[2])
    ax.plot(t_plot, np.array(df.ca),  label=labels[3], color=colors[3], linestyle=line_styles[3])
    ax.plot(t_plot, np.array(df.co),  label=labels[4], color=colors[4], linestyle=line_styles[4])

    all_y_values.extend(df.AR)
    all_y_values.extend(df.CA)
    all_y_values.extend(df.phi)
    all_y_values.extend(df.ca)
    all_y_values.extend(df.co)

    nnx = folder.split('/')[-1].replace('nnx', '')

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
plt.savefig('figs/Fig_resolution_rhytmite.png',
            format='png', dpi=300, bbox_inches='tight')

plt.close()
