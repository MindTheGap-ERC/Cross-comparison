#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import glob
import h5py

folders = [
    'data/d_phi/marlpde/Tstar_30_d_0.001',
    'data/d_phi/marlpde/Tstar_30_d_0.005',
    'data/d_phi/marlpde/Tstar_30_d_0.01',
    'data/d_phi/marlpde/Tstar_30_d_0.02',
    'data/d_phi/marlpde/Tstar_30_d_0.03',
    'data/d_phi/marlpde/Tstar_30_d_0.04',
    'data/d_phi/marlpde/Tstar_30_d_0.05',
    'data/d_phi/marlpde/Tstar_30_d_0.1',
]

labels = ["Aragonite", "Calcite", "Porosity", r"$Ca^{2+}$", r"$CO_3^{2-}$"]
line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]
colors = ['#E69F00', 'red', '#000000', '#009E73', '#0072B2']

# Depth index for time series (-1 = bottom of the system)
depth = -1

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()

Ts = 131.9 / 0.1**2

all_y_values = []

for idx, folder in enumerate(folders):
    h5_files = glob.glob(f'{folder}/*.hdf5')
    if not h5_files:
        print(f"Warning: No .hdf5 file found in {folder}")
        continue
    filename = h5_files[0]
    ax = axes[idx]

    sol = h5py.File(filename, 'r')
    df = np.array(sol['data'])
    dt = np.array(sol['times'])
    sol.close()

    t_plot = dt * Ts / 1000

    AR  = df[:, 0, depth]
    CA  = df[:, 1, depth]
    phi = df[:, 4, depth]
    Ca  = df[:, 2, depth]
    CO  = df[:, 3, depth]

    ax.plot(t_plot, AR,  label=labels[0], color=colors[0], linestyle=line_styles[0])
    ax.plot(t_plot, CA,  label=labels[1], color=colors[1], linestyle=line_styles[1])
    ax.plot(t_plot, phi, label=labels[2], color=colors[2], linestyle=line_styles[2])
    ax.plot(t_plot, Ca,  label=labels[3], color=colors[3], linestyle=line_styles[3])
    ax.plot(t_plot, CO,  label=labels[4], color=colors[4], linestyle=line_styles[4])

    all_y_values.extend(AR)
    all_y_values.extend(CA)
    all_y_values.extend(phi)
    all_y_values.extend(Ca)
    all_y_values.extend(CO)

    d_phi = folder.split('_d_')[-1]

    ax.set_title(rf'$D_\phi$ = {d_phi}', fontsize=10)
    ax.set_xlabel('Time [ky]')
    ax.set_ylabel('Concentration/Porosity')

    panel_label = chr(97 + idx)
    ax.text(0.02, 0.98, f'({panel_label})', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left')

    if idx == 7:
        ax.legend(loc='lower right', fontsize=8)

y_min = min(all_y_values)
y_max = max(all_y_values)
y_margin = (y_max - y_min) * 0.05
for ax in axes:
    ax.set_ylim(y_min - y_margin, y_max + y_margin)

plt.tight_layout()
plt.savefig('figs/Fig_A_d_phi_marlpde.png',
            format='png', dpi=300, bbox_inches='tight')

plt.close()
