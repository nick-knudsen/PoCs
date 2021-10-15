import numpy as np
import operator
import plotly.graph_objs as go
import plotly.express as px
from scipy import stats


# Get the denominator of the probability mass function
total = 0
integers = list(range(1, 10_000_000))
for k in integers:
    total += k**-2.5

def better_get_power_law_value(*args):
    gamma = 5/2
    return (1 - np.random.rand()) ** (-1 / (gamma - 1))

def get_max_vals(samples: int, sample_size: int):
    max_vals = np.zeros(samples)
    vfunc = np.vectorize(better_get_power_law_value)

    for i in range(samples):
        values = np.fromfunction(vfunc, (sample_size,), dtype=float)
        # values = np.fromfunction(lambda i: (1 - np.random.rand() ** (-1 / (5/2 - 1))), (sample_size,), dtype=float)
        max_vals[i] = np.max(values)
        if sample_size == 1000000:
            print(i)

    return max_vals

data = dict()
for exp in range(1, 7):
    data[10**exp] = list()

for sample_size in data.keys():
    print(f'getting max vals for {sample_size}')
    data[sample_size] = get_max_vals(1000, sample_size)

# # plots for each
for sample_size in data.keys():
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(1000)), y=data[sample_size], name=f"Max Values for Sample Sizes of {sample_size:,}", mode='markers'))
    fig.show()

# fit....
x = []
y = []
for sample_size in data.keys():
    mean_value = np.log10(np.mean(data[sample_size]))
    log_size = np.log10(sample_size)
    y.append(mean_value)
    x.append(log_size)

regr = stats.linregress(x, y)

x = np.asarray(x)
y = np.asarray(y)

fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, name='Points', mode='markers'))
fig.add_trace(go.Scatter(x=x, y=regr.intercept + x*regr.slope, name='Regression Fit', mode='lines'))
fig.update_layout(title=f'Scaling of Expected Largest Value as a Function of Sample Size N. Slope of {regr.slope}')
fig.show()
