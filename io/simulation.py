from typing import Union
import numpy as np

# Linear algebra types:
Vector = list[float]
Matrix = list[list[float]]

# Number of samples:
n: int = 10
# Number of loci:
p: int = 42
# Sample marker effects:
b: Vector = np.reshape(np.random.uniform(-10, 10, p), (p, 1))
# Simulate linkage disequilibrium by sample from a multivariate normal distribution:
C: Matrix = np.matmul(b, b.T)
X: Matrix = np.random.multivariate_normal(np.reshape(b, (p,)), C, n, check_valid='raise')
X.shape