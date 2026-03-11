#!/usr/bin/env python3
# Скрипт для визуализации аттрактора Лоренца (увеличенные графики)

import numpy as np
import matplotlib.pyplot as plt
import os

# Путь к данным
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

def load_data(filename):
    """Загружает CSV с колонками t,x,y,z"""
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"Файл {path} не найден")
        return None, None, None, None
    try:
        data = np.loadtxt(path, delimiter=',', skiprows=1)
        return data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    except Exception as e:
        print(f"Ошибка загрузки {filename}: {e}")
        return None, None, None, None

# Загружаю данные (например, для RK4)
t, x, y, z = load_data('rk4.csv')
if t is None:
    exit(1)

# Нормировка времени для цвета
time_norm = (t - t[0]) / (t[-1] - t[0])

# Настройки крупных графиков
plt.rcParams.update({'font.size': 14})   # увеличенный шрифт

# Создаю фигуру с тремя подграфиками, большого размера
fig, axes = plt.subplots(1, 3, figsize=(24, 7))  # ширина 24 дюйма, высота 7
fig.suptitle('Аттрактор Лоренца (σ=10, ρ=28, β=8/3)', fontsize=18)

# Проекции и подписи
projs = [('x', 'y', x, y), ('x', 'z', x, z), ('y', 'z', y, z)]

for ax, (xl, yl, u, v) in zip(axes, projs):
    sc = ax.scatter(u, v, c=time_norm, cmap='viridis', s=1.5, alpha=0.8, edgecolors='none')
    ax.set_xlabel(xl, fontsize=14)
    ax.set_ylabel(yl, fontsize=14)
    ax.set_title(f'Проекция {xl}-{yl}', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

# Цветовая шкала
cbar = fig.colorbar(sc, ax=axes[-1], label='Время')
cbar.ax.yaxis.label.set_size(14)
cbar.ax.tick_params(labelsize=12)

plt.tight_layout()

# Сохраняю с высоким разрешением
out_path = os.path.join(DATA_DIR, 'attractor_large.png')
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"График сохранён: {out_path} (размер {24*300}x{7*300} пикселей)")

plt.show()