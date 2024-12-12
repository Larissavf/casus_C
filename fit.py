import numpy as np
from matplotlib import pyplot as plt

def solve_rungekutta(t, a,y0, ymax,):
    y = y0                              # Initial condition
    steps = int(abs(t) / 0.01) + 1      # Number of ~0.01 time steps
    dt = t / steps

    for step in range(steps):
        dydt1 = a * y * np.log(ymax/y)
        y1 = y + dydt1 * 0.5 * dt
        dydt2 = a * y1 * np.log(ymax/y1)
        y2 = y + dydt2 * 0.5 * dt
        dydt3 = a * y2 * np.log(ymax/y2)
        y3 = y + dydt3 * dt
        dydt4 = a * y3 * np.log(ymax/y3) 
        dydt = (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0
        y += dydt * dt

    return y

# Original data
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

# Scale the data_ts to [-1, 1]
data_ts = np.array(data_ts)
min_val = np.min(data_ts)
max_val = np.max(data_ts)
scaled_ts = 2 * (data_ts - min_val) / (max_val - min_val) - 1

data_ys = np.array(data_ys)
scaled_ys = data_ys/np.max(data_ys)

print(scaled_ys)


print(scaled_ts)

# Plot scaled data and model predictions
ts = np.linspace(-1, 1, 600)  # Scaled time series for model
def mean_squared_error(data_ts, data_ys, **params):
    N = len(data_ts)
    total = 0.0
    for i in range(N):
        error = data_ys[i] - solve_rungekutta(data_ts[i], **params)
        total += error * error
    return total / N

params = {'a': 3.533, 'y0': 5.788, 'ymax': 7.479}
ys = [solve_rungekutta(t, **params) for t in ts]

plt.plot(ts, ys, '-', label='model')
plt.plot(scaled_ts, scaled_ys, 'o', label='data')
plt.axhline(0.0, color='k'); plt.axvline(0.0, color='k')
plt.grid(True); plt.legend()
plt.xlabel('$t$'); plt.ylabel('$y(t)$')
plt.show()

# Optimize parameters
deltas = {key: 1.0 for key in params}
params = {'a': 3.533, 'y0': 5.788, 'ymax': 7.479}
mse = mean_squared_error(scaled_ts, scaled_ys, **params)
counter = 0

while max(abs(delta) for delta in deltas.values()) > 1e-6:
    counter += 1
    print("loop:", counter)
    for key in params:
        new_params = params.copy()
        new_params[key] = params[key] + deltas[key]
      #   for param in new_params:
      #       if new_params[param] < 0:
      #           new_params[param] = 1e-8
        new_mse = mean_squared_error(scaled_ts, scaled_ys, **new_params)
        if new_mse < mse:
            params = new_params
            mse = new_mse
            deltas[key] *= 1.2
            continue

        new_params[key] = params[key] - deltas[key]
      #   for param in new_params:
      #       if new_params[param] < 0:
      #           new_params[param] = 1e-8
        new_mse = mean_squared_error(scaled_ts, scaled_ys, **new_params)
        if new_mse < mse:
            params = new_params
            mse = new_mse
            deltas[key] *= -1.0
            continue

        deltas[key] *= 0.5

print('Optimized parameters:')
for key, val in params.items():
    print(f'* {key:>2s} = {val:9.6f}')

aic = len(data_ts) * np.log(mean_squared_error(scaled_ts, scaled_ys, **params)) + 2*len(params)
print(aic)
