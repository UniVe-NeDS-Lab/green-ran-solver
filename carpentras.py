import csv
import matplotlib.pyplot as plt
import numpy as np

data = []
with open('Carpentras.csv', newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        data.append(row)

print(len(data))

def gen_hourly():
    hourly = []
    tmp = 0
    cnt = 0
    prev = 'zero'
    for d in data[2:]:
        if int(d[7]) > 0:
            if prev == 'normal':
                hourly.append(int(d[7]))
            prev = 'normal'
        elif int(d[7]) == 0:
            if prev == 'normal':
                prev = 'zero'
                hourly = hourly[:-1]

    hour_max = max(hourly)
    rain_max = hour_max/5
    cloudy_max = hour_max/2
    print(hour_max)
    print(cloudy_max)
    print(rain_max)

    sc = []
    sr = []
    cs = []
    cr = []
    rs = []
    rc = []

    hour_pass = 0
    prev = 's'
    for d in hourly:
        curr = ''
        if d < rain_max:
            curr = 'r'
        elif d < cloudy_max:
            curr = 'c'
        else:
            curr = 's'
        if prev == 's':
            if curr == 'c':
                sc.append(hour_pass)
                hour_pass = 0
            elif curr == 'r':
                sr.append(hour_pass)
                hour_pass = 0
        elif prev == 'c':
            if curr == 's':
                cs.append(hour_pass)
                hour_pass = 0
            elif curr == 'r':
                cr.append(hour_pass)
                hour_pass = 0
        elif prev == 'r':
            if curr == 's':
                rs.append(hour_pass)
                hour_pass = 0
            elif curr == 'c':
                rc.append(hour_pass)
                hour_pass = 0
        hour_pass += 1
        prev = curr

    '''print(sc)
    print(sr)
    print(cs)
    print(cr)
    print(rs)
    print(rc)'''

    print("Rate of Sunny -> Cloudy: "+str(sum(sc)/len(sc)))
    print("Rate of Sunny -> Rain: "+str(sum(sr)/len(sr)))
    print("Rate of Cloudy -> Sunny: "+str(sum(cs)/len(cs)))
    print("Rate of Cloudy -> Rain: "+str(sum(cr)/len(cr)))
    print("Rate of Rain -> Sunny: "+str(sum(rs)/len(rs)))
    print("Rate of Rain -> Cloudy: "+str(sum(rc)/len(rc)))

def gen_daily():
    spring = [3,4,5]
    summer = [6,7,8]
    autumn = [9,10,11]
    winter = [12,1,2]
    all_seasons = spring+summer+autumn+winter

    daily = []
    tmp = 0
    cnt = 0
    for d in data[2:]:
        if int(d[0].split('/')[0]) in all_seasons:
            if int(d[7]) > 0:
                tmp += int(d[7])
                cnt += 1
            if d[1] == '24:00':
                avg = tmp/cnt
                daily.append(avg)
                tmp = 0
                cnt = 0

    day_max = max(daily)
    rain_max = day_max/5
    cloudy_max = day_max/2

    sc = []
    sr = []
    cs = []
    cr = []
    rs = []
    rc = []

    day_pass = 0
    prev = 's'
    for d in daily:
        curr = ''
        if d < rain_max:
            curr = 'r'
        elif d < cloudy_max:
            curr = 'c'
        else:
            curr = 's'
        if prev == 's':
            if curr == 'c':
                sc.append(day_pass)
                day_pass = 0
            elif curr == 'r':
                sr.append(day_pass)
                day_pass = 0
        elif prev == 'c':
            if curr == 's':
                cs.append(day_pass)
                day_pass = 0
            elif curr == 'r':
                cr.append(day_pass)
                day_pass = 0
        elif prev == 'r':
            if curr == 's':
                rs.append(day_pass)
                day_pass = 0
            elif curr == 'c':
                rc.append(day_pass)
                day_pass = 0
        day_pass += 1
        prev = curr

    '''print(sc)
    print(sr)
    print(cs)
    print(cr)
    print(rs)
    print(rc)'''

    print("sc_rate = "+str(sum(sc)/len(sc)))
    print("sr_rate = "+str(sum(sr)/len(sr)))
    print("cs_rate =  "+str(sum(cs)/len(cs)))
    print("cr_rate =  "+str(sum(cr)/len(cr)))
    print("rs_rate =  "+str(sum(rs)/len(rs)))
    print("rc_rate =  "+str(sum(rc)/len(rc)))

#gen_hourly()
gen_daily()
