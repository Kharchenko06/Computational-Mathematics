import os

ROOT     = os.path.join(os.path.dirname(__file__), '..')
DATA_DIR = os.path.join(ROOT, 'data')

# Читаем параметры из data/run_params.cfg — файл создаётся программой при каждом запуске.
# Это гарантирует, что заголовки графиков всегда соответствуют реально использованным параметрам,
# независимо от того, были ли они введены вручную или загружены из lorenz.cfg.
_RUN_PARAMS = os.path.join(DATA_DIR, 'run_params.cfg')

def load_cfg():
    defaults = {
        'sigma': 10.0, 'rho': 28.0, 'beta': 8.0/3.0,
        'x0': 1.0, 'y0': 1.0, 'z0': 1.0,
        't_start': 0.0, 't_end': 50.0,
        'h_rk4': 0.01, 'h_rk38': 0.01,
        'rtol': 1e-7, 'atol': 1e-10,
    }
    path = _RUN_PARAMS
    if not os.path.exists(path):
        print(f"[config_reader] {path} не найден, используются значения по умолчанию.")
        return defaults
    with open(path) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if '=' not in line:
                continue
            key, val = line.split('=', 1)
            key, val = key.strip(), val.strip()
            if key in defaults:
                defaults[key] = float(val)
    return defaults

def make_title(cfg, prefix=''):
    beta_str = f"{cfg['beta']:.4g}"
    y0_str   = f"[{cfg['x0']:.4g}, {cfg['y0']:.4g}, {cfg['z0']:.4g}]"
    t_str    = f"[{cfg['t_start']:.4g}, {cfg['t_end']:.4g}]"
    parts = [
        f"σ={cfg['sigma']:.4g}",
        f"ρ={cfg['rho']:.4g}",
        f"β={beta_str}",
        f"y₀={y0_str}",
        f"t∈{t_str}",
    ]
    sep = '  '
    line = sep.join(parts)
    return f"{prefix}  {line}" if prefix else line

def load_trajectory(filename):
    import numpy as np
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"Файл не найден: {path}")
        return None, None, None, None
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def load_error(filename):
    import numpy as np
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"Файл не найден: {path}")
        return None, None
    data = np.loadtxt(path, delimiter=',', skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]  # t, abs_error