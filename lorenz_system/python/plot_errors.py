#!/usr/bin/env python3
# Скрипт для визуализации ошибок методов интегрирования

import numpy as np
import matplotlib.pyplot as plt
import os

# Путь к данным
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

def load_data(filename):
    """Загружает CSV с колонками t, error"""
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"Файл {path} не найден")
        return None, None
    try:
        data = np.loadtxt(path, delimiter=',', skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data[:, 0], data[:, 1]
    except Exception as e:
        print(f"Ошибка загрузки {filename}: {e}")
        return None, None

err_files = {
    'RK4':  ('error_rk4.csv',  'steelblue'),
    'RK38': ('error_rk38.csv', 'darkorange'),
    'RK23': ('error_rk23.csv', 'forestgreen'),
}

plt.rcParams.update({'font.size': 14})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9))
fig.suptitle('Ошибки методов относительно эталона (RK4, h=1e-5)', fontsize=18)

for name, (fname, col) in err_files.items():
    t_e, err = load_data(fname)
    if t_e is None:
        continue
    ax1.plot(t_e, err, color=col, lw=1.2, label=name)
    ax2.semilogy(t_e, np.where(err > 0, err, 1e-20), color=col, lw=1.2, label=name)

for ax, title in [(ax1, 'Линейная шкала'), (ax2, 'Логарифмическая шкала')]:
    ax.set_ylabel('||ошибка||', fontsize=14)
    ax.set_xlabel('t', fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=12)

plt.tight_layout()

out_path = os.path.join(DATA_DIR, 'errors.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"График сохранён: {out_path}")

plt.show()