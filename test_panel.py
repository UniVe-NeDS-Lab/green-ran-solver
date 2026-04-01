from base import multiple_bs
import numpy as np
from scipy.linalg import null_space
import matplotlib.pyplot as plt
import time
from scipy.linalg import svd
from scipy.sparse import lil_matrix
import time
import csv

header = []
meas = ['Produced','Wasted','Consumed','Bought',
        'Wasted - Sunny', 'Wasted - Cloudy', 'Wasted - Rain',
        'Produced - Sunny', 'Produced - Cloudy', 'Produced - Rain',
        'Consumed - Sunny', 'Consumed - Cloudy', 'Consumed - Rain',
        'Bought - Sunny', 'Bought - Cloudy', 'Bought - Rain',
        'Battery Empty','Battery Full',
        'Empty - Sunny', 'Empty - Cloudy', 'Empty - Rain',
        'Full - Sunny', 'Full - Cloudy', 'Full - Rain',
        'Sunny','Cloudy','Rain']

tc = 6
panels = [(tc/4)*3,tc/4]

tb = 100
batteries = [75,25]

for bid in range(len(batteries)):
    for m in meas:
        header.append("Battery "+str(bid)+" "+m)

for m in meas:
    header.append("Total "+m)

print(header)

bs_number = 4
# energy request size
ks = [1,2,4,2]
# enery request rate
lams = [0.3,0.2,0.5, 0.2]

sc = 1/(4.147058823529412*1)
sr = 1/(3.7333333333333334*1)
cs = 1/(1.105263157894737*1)
cr = 1/(1.2*1)
rs = 1/(1.2692307692307692*1)
rc = 1/(1.6428571428571428*1)
wr = np.array([[0, sc, sr],[cs, 0, cr], [rs, rc, 0]])

panel_map = [[0],[1]] #panel to battery [[0,1],[2]]

pb_dist = [0,0] #distance panel -> battery

bs_map = [[0,1,2],[3]] #battery to bs [[0],[3],[1,2]]
bb_dist = [1.0,0.7,1.7,0]

beta1 = 0.1 #lost energy entering battery
beta2 = 0.1 #            exiting

panelss = np.linspace(4,20,21,dtype=float)
print(panelss)

csv_res = []

header = ['Total Panel Capacity'] + header

for pan in panelss:
    tc = pan
    panels = [(tc/4)*3,tc/4]
    print(panels)
    res = multiple_bs(tc,panels,tb,batteries,bs_number,ks,lams,wr,panel_map,pb_dist,bs_map,bb_dist,beta1,beta2)
    csv_res.append([pan]+res)

with open("Test_Case2_panel.csv", "w", newline="") as f:
    writer = csv.writer(f)
    # Optional header
    writer.writerow(header)
    # Write rows
    for i in range(len(csv_res)):
        writer.writerow(csv_res[i])