from typing import Union
import numpy as np
import pandas as pd
import scipy as sp

# Linear algebra types:
Vector = list[float]
Matrix = list[list[float]]

# Number of samples:
n: int = 1_00
# Number of loci:
p: int = 42_000
# Number of equally size chromosomes
m: int = 7
# Sample marker effects:
b: Vector = np.zeros(p)
b: Vector = np.reshape(np.random.uniform(-10, 10, p), (p, 1))
# Simulate linkage disequilibrium by sample from a multivariate normal distribution:
C: Matrix = np.matmul(b, b.T)
for j in range(p):
    C[j][j] += 1.00

# Numpy way
X: Matrix = np.random.multivariate_normal(np.reshape(b, (p,)), C, n, check_valid='raise')
X.shape



# Scipy way
dist = sp.stats.multivariate_normal(mean=np.reshape(b, (p,)), cov=C)

X = np.empty((n, p))
for i in range(n):
    X[i, :] = dist.rvs(size=1)

X