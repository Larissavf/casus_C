def solve_rungekutta(t, a, b, y0, ymax):
    y = y0                              # Beginconditie
    steps = int(abs(t) / 0.01) + 1      # Hoeveel tijdstappen van ~0.01 zijn nodig
    dt = t / steps
    for step in range(steps):
        # Tijdelijke stappen:
        dydt1 = a * y * (ymax**b - y**b)            # Differentiaalvergelijking stap 1
        y1 = y + dydt1 * 0.5 * dt
        dydt2 = a * y1 * (ymax**b - y1**b)             # Differentiaalvergelijking stap 2
        y2 = y + dydt2 * 0.5 * dt
        dydt3 = a * y2 * (ymax**b - y2**b)                  # Differentiaalvergelijking stap 3
        y3 = y + dydt3 * dt
        dydt4 = a * y3 * (ymax**b - y3**b)                  # Differentiaalvergelijking stap 4
        # Definitieve stap:
        dydt = (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0
        y += dydt * dt
    return y

params = {
    'a':  0.859,
    'b':  0.121,
    'y0': 0.000755,
    "ymax": 7.45
}
ts = [i / 10 for i in range(600)]
ys = [solve_rungekutta(t, **params) for t in ts]

from matplotlib import pyplot as plt

data_ts = [
     3.46,  4.58,  5.67,  6.64,  7.63,  8.41,  9.32, 10.27, 11.19,
    12.39, 13.42, 15.19, 16.24, 17.23, 18.18, 19.29, 21.23, 21.99,
    24.33, 25.58, 26.43, 27.44, 28.43, 30.49, 31.34, 32.34, 33.00,
    35.20, 36.34, 37.29, 38.50, 39.67, 41.37, 42.58, 45.39, 46.38,
    48.29, 49.24, 50.19, 51.14, 52.10, 54.00, 56.33, 57.33, 59.38,
]
data_ys = [
    0.0158, 0.0264, 0.0326, 0.0445, 0.0646, 0.0933, 0.1454, 0.2183, 0.2842,
    0.4977, 0.6033, 0.8441, 1.2163, 1.4470, 2.3298, 2.5342, 3.0064, 3.4044,
    3.2046, 4.5241, 4.3459, 5.1374, 5.5376, 4.8946, 5.0660, 6.1494, 6.8548,
    5.9668, 6.6945, 6.6395, 6.8971, 7.2966, 7.2268, 6.8815, 8.0993, 7.2112,
    7.0694, 7.4971, 6.9974, 6.7219, 7.0523, 7.1095, 7.0694, 8.0562, 7.2268, 
]

plt.plot(ts, ys, '-', label='model')
plt.plot(data_ts, data_ys, 'o', label='data')
plt.axhline(0.0, color='k'); plt.axvline(0.0, color='k')
plt.grid(True); plt.legend()
plt.xlabel('$t$'); plt.ylabel('$y(t)$')
plt.show()

def mean_squared_error(data_ts, data_ys, **params):
    N = len(data_ts)
    total = 0.0
    for i in range(N):
        error = data_ys[i] - solve_rungekutta(data_ts[i], **params)
        total += error * error
    return total / N

print(mean_squared_error(data_ts, data_ys, **params))

# Initialisatie
params = {
    'a':  0.859,
    'b':  0.121,
    'y0': 0.000755,
    'ymax': 7.45
}
deltas = {key: 1.0 for key in params}

# Herhaaldelijke aanpassing
mse = mean_squared_error(data_ts, data_ys, **params)
counter = 0
while max(abs(delta) for delta in deltas.values()) > 0.001:
    counter += 1
    print("loop:", counter)
    for key in params:
        new_params = params.copy()
        # Probeer de betreffende parameter the verhogen
        new_params[key] = params[key] + deltas[key]
        for param in new_params:
            if new_params[param] < 0:
                new_params[param] = 0
        new_mse = mean_squared_error(data_ts, data_ys, **new_params)
        if new_mse < mse:
            params = new_params
            mse = new_mse
            deltas[key] *= 1.2
            continue
        # Probeer de betreffende parameter the verlagen
        new_params[key] = params[key] - deltas[key]
        for param in new_params:
            if new_params[param] < 0:
                new_params[param] = 0
        new_mse = mean_squared_error(data_ts, data_ys, **new_params)
        if new_mse < mse:
            params = new_params
            mse = new_mse
            deltas[key] *= -1.0
            continue
        # Verklein de stapgrootte
        deltas[key] *= 0.5

print('Optimale parameters:')
for key, val in params.items():
    print(f'* {key:>2s} = {val:9.6f}')
