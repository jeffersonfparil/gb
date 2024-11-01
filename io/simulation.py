from typing import Union
import numpy as np
import polars as pl
import scipy as sp

# Type aliases
VectorF = list[float]
MatrixF = list[list[float]]
VectorI = list[int]
MatrixI = list[list[int]]
VectorB = list[bool]
MatrixB = list[list[bool]]
VectorS = list[str]
MatrixS = list[list[str]]
# Number of genotypes:
n: int = 1_00
if not (n >= 1):
    raise Exception('The number of genotypes (n) should be non-zero.')

# Number of loci:
p: int = 42_693
if not (p >= 1):
    raise Exception('The number of loci (p) should be non-zero.')

# Number of equally size chromosomes
m: int = 5
if not (m >= 1) and (m <= p):
    raise Exception('The number of equally sized chromosomes (m) should range from one to the total number of loci.')

# Fraction of loci with effects
f = 0.01
if not (f >= 0) and (f <= 1):
    raise Exception('The fraction of loci with effects (f) should range from zero to one.')

# Number of loci with effects
a: int = round(f * p)
if not a > 0:
    raise Exception('No markers with effects. Please modify p (total number of loci) and/or f (fraction of loci with effects).')

# Constants
B_MU: float = 10.00 ### mean of the normally distributed non-zero loci effects
B_SD: float = 2.00 ### standard deviation of the normally distributed non-zero loci effects
P: int = 2_987_654_321 ### maximum length of the genome in base pairs

# Sample marker effects:
b: VectorF = np.zeros(p)
b_indexes_with_effects: VectorI = np.random.randint(0, p, a)
b_effects: VectorF = np.random.normal(B_MU, B_SD, a)
b[b_indexes_with_effects] = b_effects

# Define the structure of the genome
# Number of loci per chromosome
chr_nloci_base: int = np.round(p / m)
chr_nloci: VectorI = [chr_nloci_base for i in range(m)]
if (np.sum(chr_nloci) < p):
    chr_nloci[-1] += p-np.sum(chr_nloci)
elif (np.sum(chr_nloci) > p):
    chr_nloci[-1] -= np.sum(chr_nloci)-p

# Lengths of each chromosome
chr_length_base: int = np.round(P / m)
chr_lengths: VectorI = [chr_length_base for i in range(m)]
if (np.sum(chr_lengths) < P):
    chr_lengths[-1] += P-np.sum(chr_lengths)
elif (np.sum(chr_lengths) > P):
    chr_lengths[-1] -= np.sum(chr_lengths)-P

chromosomes: VectorS = []
positions: VectorI = []
for i in range(m):
    n_loci = int(chr_nloci[i])
    chromosomes.extend(["chr_" + str(i) for j in range(n_loci)])
    tmp_positions = np.random.randint(0, int(chr_lengths[i]), n_loci)
    # Remove duplicates byt incrementing by 1 which works because np.random.randint does not include the maximum value in the choices
    tmp_positions.sort()
    for j in range(n_loci-1):
        if tmp_positions[j] == tmp_positions[j+1]:
            tmp_positions[j+1] += 1;
    positions.extend(tmp_positions)

len(chromosomes) == p
len(positions) == p



# Simulate linkage disequilibrium by sample from a multivariate normal distribution:
C: MatrixF = np.matmul(b, b.T)
for j in range(p):
    C[j][j] += 1.00

# Numpy way
X: MatrixF = np.random.multivariate_normal(np.reshape(b, (p,)), C, n, check_valid='raise')
X.shape



# Scipy way
dist = sp.stats.multivariate_normal(mean=np.reshape(b, (p,)), cov=C)

X = np.empty((n, p))
for i in range(n):
    X[i, :] = dist.rvs(size=1)

X