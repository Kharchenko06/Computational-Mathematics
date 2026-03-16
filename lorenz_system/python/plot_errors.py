#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from config_reader import load_cfg, make_title, load_error, DATA_DIR
import os

cfg = load_cfg()

err_files = {
    'RK4':  ('error_rk4.csv',  'steelblue'),
    'RK38': ('error_rk38.csv', 'darkorange'),
    'RK23': ('error_rk23.csv', 'forestgreen'),
}

plt.rcParams.update({'font.size': 14})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9))
fig.suptitle(make_title(cfg, 'Абсолютная ошибка (относительно эталона RK4, h=1e-5)'), fontsize=14)

for name, (fname, col) in err_files.items():
    t_e, err = load_error(fname)
    if t_e is None:
        continue
    ax1.plot(t_e, err, color=col, lw=1.2, label=name)
    ax2.semilogy(t_e, np.where(err > 0, err, 1e-20), color=col, lw=1.2, label=name)

for ax, title in [(ax1, 'Линейная шкала'), (ax2, 'Логарифмическая шкала')]:
    ax.set_ylabel('||ошибка||', fontsize=14)
    ax.set_xlabel('t', fontsize=14)
    ax.set_title(title, fontsize=15)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=12)

plt.tight_layout()
out = os.path.join(DATA_DIR, 'errors.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Сохранено: {out}")
plt.show()