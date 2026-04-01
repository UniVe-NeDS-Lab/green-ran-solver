from base import compute_energy
import numpy as np
from scipy.linalg import null_space
import matplotlib.pyplot as plt
import time
from scipy.linalg import svd
from scipy.sparse import lil_matrix
import time
import csv

# Parameters
M_size = 100
M = M_size-1   # number of levels
b = 3          # block size

# sun rate
sigma_s = 2
sigma_c = sigma_s/2
sigma_r = sigma_s/5
sigma = np.array([[sigma_s, 0, 0],[0, sigma_c, 0], [0, 0, sigma_r]])

#wr = np.array([[0, 0.15, 0],[0.25, 0, 0.5], [1, 0.75, 0]])
SC = 1/259200
oth = 1/86400
wr = np.array([[0, SC, 0],[oth, 0, oth], [oth, oth, 0]])

# energy request size
ks = [1,2,3,4]
# enery request rate
lams = [0.2,0.1,0.15, 0.1]

compute_energy(M_size, sigma, wr, lams, ks, True)