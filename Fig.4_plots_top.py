import matplotlib.pyplot as plt
import pandas as pd
from comparePlot import floatMarl, plotTemporalRhy


# Create figure with 2x2 subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

############# FORTRAN RUNS ##############

# Common settings
depth_ind = 1
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
ax2.legend(loc='lower left')
ax2.text(0.02, 0.98, '(b)', transform=ax2.transAxes,
        fontsize=12, fontweight='bold', va='top')

############# RHYTHMITE RUNS ##############

# Plot third panel (bottom left - Fig 4c)
plotTemporalRhy("data/replication/rhythmite/Fig._4a/solution_x_000050.ascii", ax=ax3)
ax3.set_xlabel('Time [ky]')
ax3.set_ylabel('Concentration/Porosity')
ax3.text(0.02, 0.98, '(c)', transform=ax3.transAxes,
        fontsize=12, fontweight='bold', va='top')

# Plot fourth panel (bottom right - Fig 4d)
plotTemporalRhy("data/replication/rhythmite/Fig._4b/solution_x_000050.ascii", ax=ax4)
ax4.set_xlabel('Time [ky]')
ax4.set_ylabel('Concentration/Porosity')
ax4.legend(loc='lower left')
ax4.text(0.02, 0.98, '(d)', transform=ax4.transAxes,
        fontsize=12, fontweight='bold', va='top')

# Adjust layout and save
plt.tight_layout()
plt.savefig('Fig.4_replication_top.png', format='png', bbox_inches='tight')
plt.close()

