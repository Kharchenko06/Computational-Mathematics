#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from config_reader import load_cfg, make_title, load_trajectory, DATA_DIR
import os

cfg = load_cfg()
t, x, y, z = load_trajectory('rk4.csv')
if t is None:
    exit(1)

time_norm = (t - t[0]) / (t[-1] - t[0])

plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
fig.suptitle(make_title(cfg, 'Аттрактор Лоренца'), fontsize=16)

for ax, (xl, yl, u, v) in zip(axes, [('x','y',x,y), ('x','z',x,z), ('y','z',y,z)]):
    sc = ax.scatter(u, v, c=time_norm, cmap='viridis', s=1.5, alpha=0.8, edgecolors='none')
    ax.set_xlabel(xl, fontsize=14)
    ax.set_ylabel(yl, fontsize=14)
    ax.set_title(f'Проекция {xl}-{yl}', fontsize=15)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

cbar = fig.colorbar(sc, ax=axes[-1], label='Время')
cbar.ax.yaxis.label.set_size(13)

plt.tight_layout()
out = os.path.join(DATA_DIR, 'attractor_large.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Сохранено: {out}")
plt.show()