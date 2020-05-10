import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def p2rgb(t):
    def f(theta):
        return 0.45*(1+np.cos(theta))
    return (f(t), f(t - 2*np.pi/3), f(t + 2*np.pi/3))

data = pd.read_csv("init.csv")
plt.scatter(data.iloc[:, 0].values, data.iloc[:, 1].values, c=[p2rgb(t) for t in data.iloc[:, 2].values])
plt.title('Initial positions and phases')
plt.savefig('init.png')
plt.show()

data = pd.read_csv("final.csv")
plt.scatter(data.iloc[:, 0].values, data.iloc[:, 1].values, c=[p2rgb(t) for t in data.iloc[:, 2].values])
plt.title('Final positions and phases')
plt.savefig('final.png')
plt.show()
