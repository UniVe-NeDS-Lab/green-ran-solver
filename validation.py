from base import create_block_qbd, create_block_lite, linear_reverse
import numpy as np
from scipy.linalg import null_space
import matplotlib.pyplot as plt
import time
from scipy.linalg import svd
from scipy.sparse import lil_matrix
import time
import csv

# Parameters
M = 100       # number of levels
M_size = M+1
b = 3       # block size

# enery request rate
lams = [0.2,0.1,0.15, 0.1]

# sun rate
sigma_s = 1
sigma_c = 0.2
sigma_r = 0.1
sigma = np.array([[sigma_s, 0, 0],[0, sigma_c, 0], [0, 0, sigma_r]])

# energy request size
ks = [1,2,3,4]

#wr = np.array([[0, 0.15, 0],[0.25, 0, 0.5], [1, 0.75, 0]])
SC = 1/259200
oth = 1/86400
wr = np.array([[0, SC, 0],[oth, 0, oth], [oth, oth, 0]])

start = time.perf_counter()
# Create block QBD matrix
Q = create_block_qbd(M_size, b, sigma, lams, ks, wr)

a_mod = Q.toarray()
a_mod = np.transpose(a_mod)
a_mod[1,:] = np.ones(a_mod.shape[0])
b_vector = np.zeros(a_mod.shape[0])
b_vector[1] = 1
b_vector
raw_pi = np.linalg.solve(a_mod,b_vector)
trad_pis = []
for i in range(M_size):
    trad_pis.append(raw_pi[b*i:b*(i+1)])

end = time.perf_counter()
print(f"Elapsed time (traditional solver): {end - start:.6f} seconds")

start = time.perf_counter()

# Create block QBD matrix
Q = create_block_lite(M_size, b, sigma, lams, ks, wr)
end = time.perf_counter()
#print(f"Elapsed time blocks init: {end - start:.6f} seconds")

max_k = max(ks)
P00 = Q[:b*max_k, :b*max_k].toarray()
A0 = Q[:b*max_k, b*max_k:2*b*max_k].toarray()
A1 = Q[b*max_k:2*b*max_k,b*max_k:2*b*max_k].toarray()
A2 = Q[b*max_k:2*b*max_k,:b*max_k].toarray()
lb = 0
if ((M+1)%max_k) == 0:
    lb = int(((M+1)/max_k)-1)
else:
    lb = ((M+1)//max_k)
PMM = Q[b*max_k*lb:,b*max_k*lb:].toarray()
P01 = Q[b*max_k*(lb-1):b*max_k*lb,b*max_k*lb:].toarray()
P10 = Q[b*max_k*lb:,b*max_k*(lb-1):b*max_k*lb].toarray()
if ((M+1)//max_k) < 2:
    A1 = P00
    A2 = P10

new_M = 0
if ((M+1)%max_k) == 0:
    new_M = int(((M+1)/max_k)-1)
else:
    new_M = ((M+1)//max_k)

end = time.perf_counter()
#print(f"Elapsed time blocks generating: {end - start:.6f} seconds")
pis = linear_reverse(new_M,b,P00,P01,A0,A1,A2,PMM,P10)

end = time.perf_counter()
print(f"Elapsed time (fast solver): {end - start:.6f} seconds")

for i, (p1, p2) in enumerate(zip(trad_pis, pis)):
    if np.allclose(p1, p2, atol=1e-16):
        print(f"Level {i}: ✅ match")
    else:
        print(f"Level {i}: ❌ mismatch")
        print(f"  Normal: {p1}")
        print(f"  Folded: {p2}")
pass