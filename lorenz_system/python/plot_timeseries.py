#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from config_reader import load_cfg, make_title, load_trajectory, DATA_DIR
import os

cfg = load_cfg()

methods = [
    (f"RK4 (h={cfg['h_rk4']})",   'rk4.csv',  'steelblue',  '-'),
    (f"RK3/8 (h={cfg['h_rk38']})", 'rk38.csv', 'darkorange', '--'),
    ('RK23 (адаптивный)',           'rk23.csv', 'forestgreen',':'),
]

T_LIMIT = cfg['t_end']

plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(3, 1, figsize=(16, 9), sharex=True)
fig.suptitle(make_title(cfg, 'Система Лоренца: x(t), y(t), z(t)'), fontsize=15)

for label, fname, color, ls in methods:
    t, x, y, z = load_trajectory(fname)
    if t is None:
        continue
    mask = t <= T_LIMIT
    for ax, comp in zip(axes, [x[mask], y[mask], z[mask]]):
        ax.plot(t[mask], comp, color=color, linestyle=ls, lw=1.2, label=label, alpha=0.85)

for ax, lbl in zip(axes, ['x(t)', 'y(t)', 'z(t)']):
    ax.set_ylabel(lbl, fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

axes[0].legend(fontsize=12, loc='upper right')
axes[-1].set_xlabel('t', fontsize=14)

plt.tight_layout()
out = os.path.join(DATA_DIR, 'timeseries.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Сохранено: {out}")
plt.show()