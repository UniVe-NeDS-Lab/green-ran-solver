import numpy as np
from scipy.linalg import null_space
import matplotlib.pyplot as plt
import time
from scipy.linalg import svd
from scipy.sparse import lil_matrix
import time
import csv

def create_block_qbd(M, b, sigma, lams, ks, weather_rate):
    total_size = M * b
    Q = lil_matrix((total_size, total_size))

    # Define block matrices
    P = weather_rate
    A1 = sigma          # arrival: go up a level
    Ad = P              # stay at same level
    B0 = Ad.copy()                # top level can only go up or stay
    Bn = Ad.copy()               # bottom level can only go down

    for i in range(M):
        row = i * b
        row_end = row + b

        # Diagonal blocks
        col = i * b
        col_end = col + b
        if i == 0:
            Q[row:row_end, col:col_end] = B0
        elif i == M - 1:
            Q[row:row_end, col:col_end] = Bn
        else:
            Q[row:row_end, col:col_end] = Ad

        # sub diagonal
        for j in range (len(lams)):
            if i >= ks[j]:
                col = (i - ks[j]) * b
                col_end = col + b
                Q[row:row_end, col:col_end] += lams[j]*np.eye(b)
            elif i > 0 :
                col = 0
                col_end = col + b
                Q[row:row_end, col:col_end] += lams[j]*np.eye(b)

        # Super-diagonal
        if i < M - 1:
            col = (i + 1) * b
            col_end = col + b
            Q[row:row_end, col:col_end] = A1

    # Ensure that the row sums are 0 for an irreducible matrix
    for i in range(total_size):
        Q[i, i] -= np.sum(Q[i, :])

    return Q

def create_block_lite(M, b, sigma, lams, ks, weather_rate):
    total_size = M * b
    Q = lil_matrix((total_size, total_size))
    max_k = max(ks)

    # Define block matrices
    P = weather_rate
    A1 = sigma          # arrival: go up a level
    Ad = P              # stay at same level
    B0 = Ad.copy()                # top level can only go up or stay
    Bn = Ad.copy()               # bottom level can only go down

    #for i in range(M):
    levels = list(set(list(range(2*max_k))+list(range(M-(2*max_k),M))))
    for i in levels:
        if i >= M or i < 0:
            continue
        row = i * b
        row_end = row + b
        row_end_new = row_end + b

        # Diagonal blocks
        col = i * b
        col_end = col + b
        col_end_new = col_end + b
        if i == 0:
            Q[row:row_end, col:col_end] = B0
        elif i == M - 1:
            Q[row:row_end, col:col_end] = Bn
        else:
            Q[row:row_end, col:col_end] = Ad

        # sub diagonal
        for j in range (len(lams)):
            if i >= ks[j]:
                col = (i - ks[j]) * b
                col_end = col + b
                Q[row:row_end, col:col_end] += lams[j]*np.eye(b)
            elif i > 0 :
                col = 0
                col_end = col + b
                Q[row:row_end, col:col_end] += lams[j]*np.eye(b)

        # Super-diagonal (A1)
        if i < M - 1:
            col = (i + 1) * b
            col_end = col + b
            Q[row:row_end, col:col_end] = A1

    # Ensure that the row sums are 0 for an irreducible matrix
    for i in range(total_size):
        Q[i, i] -= np.sum(Q[i, :])

    #dense = Q.toarray()
    #np.savetxt("q.csv", dense, delimiter=",")

    return Q

def linear_reverse(M,b,P00,P01,A0,A1,A2,PMM,P10):
    Us = [None] * (M+1)
    pis = [None] * (M+1)
    Rs = [None] * (M+1)

    Us[M] = PMM
    Rs[M] = -(np.matmul(P01, np.linalg.inv(Us[M])))
    Us[M-1] = np.subtract(A1, np.matmul(np.matmul(P01, np.linalg.inv(Us[M])), P10))
    if ((M-1) > 0):
        Rs[M-1] = -(np.matmul(A0, np.linalg.inv(Us[M-1])))
    for i in range(M-2, 0, -1):
        Us[i] = np.subtract(A1, np.matmul(np.matmul(A0, np.linalg.inv(Us[i+1])), A2))
        Rs[i] = -(np.matmul(A0, np.linalg.inv(Us[i])))
    Us[0] = np.subtract(P00, np.matmul(np.matmul(A0, np.linalg.inv(Us[1])), A2))

    pis[0] = solve_system(Us[0])
                          
    for i in range(1,M+1):
        pis[i] = np.matmul(pis[i-1], Rs[i])

    flat_pis = np.concatenate(pis)

    # Step 2: Normalize
    normalized_flat = flat_pis / np.sum(flat_pis)
    
    # Step 3: Reshape back into list of arrays
    normalized_pis = [normalized_flat[i:i+b] for i in range(0, len(normalized_flat), b)]
    return normalized_pis


def solve_system(U_m):
    a = U_m.copy()
    a = np.transpose(a)
    b_vector = np.zeros(a.shape[0])
    a[1, :] = np.ones(a.shape[0])
    b_vector[1] = 1
    raw_pi = np.linalg.solve(a,b_vector)
    raw_pi = raw_pi / np.sum(raw_pi)
    
    return raw_pi

def pad_zeros(A,target): #target = max_k * b
    n, m = A.shape
    padded = np.zeros((target, target), dtype=A.dtype)
    padded[:n, :m] = A
    return padded

def compute_energy(battery_size, sun_rate, weather_rate, lams, ks, print_res=False):
    b = 3
    M_size = battery_size
    M = M_size-1       # number of levels

    Q = create_block_lite(M_size, b, sun_rate, lams, ks, weather_rate)

    start = time.perf_counter()
    
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
    pis = linear_reverse(new_M,b,P00,P01,A0,A1,A2,PMM,P10)
    
    end = time.perf_counter()

    sigma_s = sun_rate[0][0]
    sigma_c = sun_rate[1][1]
    sigma_r = sun_rate[2][2]

    Ew = (pis[-1][0]*sigma_s) + (pis[-1][1]*sigma_c) + (pis[-1][2]*sigma_r)
    Ew_s = pis[-1][0]*sigma_s
    Ew_c = pis[-1][1]*sigma_c
    Ew_r = pis[-1][2]*sigma_r
    Ews = np.array([Ew_s,Ew_c,Ew_r])
    Ec = 0
    for i in range(len(lams)):
        Ec += (ks[i]*lams[i])

    stacked = np.vstack(pis)
    sums = np.sum(stacked, axis=0)
    Ep = (sums[0]*sigma_s) + (sums[1]*sigma_c) + (sums[2]*sigma_r)
    Ep_s = sums[0]*sigma_s
    Ep_c = sums[1]*sigma_c
    Ep_r = sums[2]*sigma_r
    Eps = np.array([Ep_s,Ep_c,Ep_r])
    Eb = Ec-(Ep-Ew)
    
    battery_empty = np.sum(pis[0])
    empty_s = pis[0][0]
    empty_c = pis[0][1]
    empty_r = pis[0][2]
    emptys = np.array([empty_s,empty_c,empty_r])
    battery_full = np.sum(pis[-1])
    full_s = pis[-1][0]
    full_c = pis[-1][1]
    full_r = pis[-1][2]
    fulls = np.array([full_s,full_c,full_r])

    if print_res:
        print('Produced: '+str(Ep))
        print('Wasted: '+str(Ew))
        print(f"Percentage Wasted: {(Ew/Ep) * 100:.2f}%")
        print('----------------')
        print('Consumed: '+str(Ec))
        print('From Battery: '+str(Ec-Eb))
        print('Bought: '+str(Eb))
        print(f"Percentage Bought: {(Eb/Ec) * 100:.2f}%")
        
        print('----------------')
        print(f"Time in sunny weather: {(sums[0]) * 100:.2f}%")
        print(f"Time in cloudy weather: {(sums[1]) * 100:.2f}%")
        print(f"Time in rainy weather: {(sums[2]) * 100:.2f}%")
        print(f"Time with empty battery: {(battery_empty) * 100:.2f}%")
        print(f"Time in full battery: {(battery_full) * 100:.2f}%")
        print('----------------')

    return (Ep, Ew, Ec, Eb, Ews, Eps, battery_empty, battery_full, emptys, fulls, sums)

def multiple_bs(tc,panels,tb,batteries,bs_number,ks,lams,wr,panel_map,pb_dist,bs_map,bb_dist,beta1,beta2):
    csv_res_tmp = []
    total_capacity = tc
    panels = panels
    
    total_battery = tb
    batteries = batteries
    
    bs_number = bs_number
    # energy request size
    ks = ks
    # enery request rate
    lams = lams
    
    wr = wr
    
    #energy_mapping = {0: {0:[(0,2.5)], 1:[(3,1.5)]}, 1: {2:[(1,0.5),(2,1)]}}
    #energy_mapping[1][2]
    
    panel_mapping = panel_map #panel to battery [[0,1],[2]]
    
    #battery to panel
    battery_number = len(batteries)
    battery_mapping = [None] * battery_number
    i = 0
    for panel in panel_mapping:
        for bid in panel:
            battery_mapping[bid] = i
        i += 1
    #print(battery_mapping)
    #battery_mapping = [0,0,1]
    
    panel_bat_distances = pb_dist #distance panel -> battery
    
    bs_mapping = bs_map #battery to bs [[0],[3],[1,2]]
    bat_bs_distances = bb_dist
    
    beta1 = beta1 #lost energy entering battery
    beta2 = beta2 #            exiting
    
    tot_res = [0] * 11
    
    for bid in range(len(battery_mapping)):
        #print("Battery "+str(bid))
        battery_size = batteries[bid]
        #print(battery_size)
        pid = battery_mapping[bid]
        panel_rate_raw = panels[pid]
        panel_rate = panel_rate_raw/len(panel_mapping[pid])
        #print(panel_rate)
        #
        sigma_s = panel_rate
        sigma_c = sigma_s/2
        sigma_r = sigma_s/5
        #lost of transportation from panel to battery
        alpha = (panel_bat_distances[bid]/1)*0.05
        sigma_s = sigma_s*(1-alpha) #transportation cost
        sigma_c = sigma_c*(1-alpha)
        sigma_r = sigma_r*(1-alpha)
        #lost when entering battery
        sigma_s = sigma_s*(1-beta1) #battery enter cost
        sigma_c = sigma_c*(1-beta1)
        sigma_r = sigma_r*(1-beta1)
        sun_rate = np.array([[sigma_s, 0, 0],[0, sigma_c, 0], [0, 0, sigma_r]])
        #print(sun_rate)
        #
        adj_lams = []
        adj_ks = []
        for bsid in bs_mapping[bid]:
            lam = lams[bsid]
            #print(lam)
            lam = (1/(1-beta2))*lam #lost when exiting battery
            alpha = (bat_bs_distances[bsid]/1)*0.05
            lam = (1/(1-alpha))*lam #lost transporting from battery to bs
            adj_lams.append(lam)
            adj_ks.append(ks[bsid])
        #print(adj_lams)
        #print(adj_ks)
        #print('-----------')
    
        res = compute_energy(battery_size, sun_rate, wr, adj_lams, adj_ks)
        csv_res_tmp += res[:4]
        csv_res_tmp += list(res[4])
        csv_res_tmp += list(res[5])
        for s in res[10]:
            csv_res_tmp += [res[2]*s]
        i = 0
        for s in res[10]:
            csv_res_tmp += [(res[2]*s)-(res[5][i]-res[4][i])]
            i+=1
        csv_res_tmp += res[6:8]
        csv_res_tmp += list(res[8])
        csv_res_tmp += list(res[9])
        csv_res_tmp += list(res[10])
        i = 0
        for r in res[:4]:
            tot_res[i] += res[i]
            i+=1
        for r in res[4:6]:
            tot_res[i] += res[i]
            i+=1
        for r in res[6:]:
            #print(res[i])
            tot_res[i] += (res[i]*(battery_size/total_battery))
            #print(tot_res[i])
            i+=1
    
    #print(tot_res)
    #print('------------')
    
    # (Ep, Ew, Ec, Eb, battery_empty, battery_full)
    Ep = tot_res[0]
    Ew = tot_res[1]
    Ec = tot_res[2]
    Eb = tot_res[3]
    battery_empty = tot_res[6]
    battery_full = tot_res[7]
    sums = tot_res[10]
    #sums = sums/battery_number
    print("Total")
    csv_res_tmp += tot_res[:4]
    csv_res_tmp += list(tot_res[4])
    csv_res_tmp += list(tot_res[5])
    for s in sums:
        csv_res_tmp += [Ec*s]
    i = 0
    for s in sums:
        csv_res_tmp += [(Ec*s)-(tot_res[5][i]-tot_res[4][i])]
        i+=1
    csv_res_tmp += tot_res[6:8]
    csv_res_tmp += list(tot_res[8])
    csv_res_tmp += list(tot_res[9])
    csv_res_tmp += list(tot_res[10])
    
    print('Produced: '+str(Ep))
    print('Wasted: '+str(Ew))
    print(f"Percentage Wasted: {(Ew/Ep) * 100:.2f}%")
    print('----------------')
    print('Consumed: '+str(Ec))
    print('From Battery: '+str(Ec-Eb))
    print('Bought: '+str(Eb))
    print(f"Percentage Bought: {(Eb/Ec) * 100:.2f}%")
    
    print('----------------')
    print(f"Time in sunny weather: {(sums[0]) * 100:.2f}%")
    print(f"Time in cloudy weather: {(sums[1]) * 100:.2f}%")
    print(f"Time in rainy weather: {(sums[2]) * 100:.2f}%")
    print(f"Time with empty battery: {(battery_empty) * 100:.2f}%")
    print(f"Time in full battery: {(battery_full) * 100:.2f}%")
    print('----------------')
    #print(csv_res_tmp)
    return csv_res_tmp