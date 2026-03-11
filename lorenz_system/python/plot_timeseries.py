#!/usr/bin/env python3
# Скрипт для построения временных рядов x(t), y(t), z(t) системы Лоренца
# Сравниваю три метода: RK4, RK3/8 и адаптивный RK23

import numpy as np
import matplotlib.pyplot as plt
import os

# Путь к данным (относительно расположения этого скрипта)
DATA_PATH = os.path.join(os.path.dirname(__file__), '..', 'data')

def load_csv(filename):
    """Загружает файл с колонками t, x, y, z (первая строка — заголовок)"""
    fullpath = os.path.join(DATA_PATH, filename)
    if not os.path.isfile(fullpath):
        print(f" [предупреждение] файл {filename} не найден, пропускаем.")
        return None, None, None, None
    try:
        # пропускаю заголовок, читаю все числа
        arr = np.loadtxt(fullpath, delimiter=',', skiprows=1)
        if arr.ndim == 1:  # если всего одна строка
            arr = arr.reshape(1, -1)
        return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    except Exception as e:
        print(f" !!! ошибка при чтении {filename}: {e}")
        return None, None, None, None

# Методы, которые сравниваем (имя, имя файла, цвет, стиль линии)
# цвета подобрал, чтобы на графике не сливались
methods = [
    ('RK4 (фикс. h=0.01)',  'rk4.csv',  'steelblue',   '-'),
    ('RK3/8 (фикс. h=0.01)', 'rk38.csv', 'darkorange',  '--'),
    ('RK23 (адаптивный)',    'rk23.csv', 'forestgreen', ':'),
]

# До какого времени показываем (чтобы хаос не забил график)
T_LIMIT = 40.0

# Настройки шрифтов
plt.rcParams.update({'font.size': 14})

# Создаю три подграфика друг под другом
fig, axes = plt.subplots(3, 1, figsize=(16, 9), sharex=True)
fig.suptitle('Система Лоренца: x(t), y(t), z(t) (σ=10, ρ=28, β=8/3)', fontsize=18)

# Цикл по методам
for label, fname, color, ls in methods:
    t, x, y, z = load_csv(fname)
    if t is None:
        continue
    # обрезаю по времени
    mask = t <= T_LIMIT
    t_ = t[mask]
    x_ = x[mask]
    y_ = y[mask]
    z_ = z[mask]

    # рисую на каждом подграфике свою компоненту
    axes[0].plot(t_, x_, color=color, linestyle=ls, linewidth=1.2, label=label, alpha=0.8)
    axes[1].plot(t_, y_, color=color, linestyle=ls, linewidth=1.2, label=label, alpha=0.8)
    axes[2].plot(t_, z_, color=color, linestyle=ls, linewidth=1.2, label=label, alpha=0.8)

# Подписи осей и прочее
axes[0].set_ylabel('x(t)', fontsize=14)
axes[1].set_ylabel('y(t)', fontsize=14)
axes[2].set_ylabel('z(t)', fontsize=14)
axes[2].set_xlabel('t', fontsize=14)

for ax in axes:
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

# Легенду размесщу вверху справа (только на первом графике, чтобы не дублировать)
axes[0].legend(fontsize=12, loc='upper right')

plt.tight_layout()

# Сохраняю в ту же папку data
out_file = os.path.join(DATA_PATH, 'timeseries.png')
plt.savefig(out_file, dpi=300, bbox_inches='tight')
print(f"График сохранён: {out_file}")

plt.show()