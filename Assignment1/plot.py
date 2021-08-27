import os
import matplotlib.pyplot as plt
import numpy as np


with open('data16.txt') as f:
    lines = (line.rstrip() for line in f) 
    z = list(line for line in lines if line)

with open('data36.txt') as f:
    lines = (line.rstrip() for line in f) 
    z = z + list(line for line in lines if line)

with open('data49.txt') as f:
    lines = (line.rstrip() for line in f) 
    z = z + list(line for line in lines if line)

with open('data64.txt') as f:
    lines = (line.rstrip() for line in f) 
    z = z + list(line for line in lines if line)

for k in range(4):
    rows, cols = (21, 5) 
    data_set = [[] for j in range(rows)]                      #   0-6 == M1   ,,,,  7-13 === M2 ,,,,  14-20 === M3

    for i in range(k*7*5*3,(k+1)*7*5*3):
        j = i%15
        D = int(i/15)%7
        if j%3 == 0 :
            data_set[D].append(z[i])
        if j%3 == 1 :
            data_set[D+7].append(z[i])
        if j%3 == 2 :
            data_set[D+14].append(z[i])


    for i in range(0,21):
        data_set[i] = np.asarray(data_set[i]).astype(float)

    ticks = ['A', 'B', 'C']

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    plt.figure()
    plt.yscale('log')
    ticks = ['16^2','32^2','64^2','128^2','256^2','512^2','1024^2']

    for i in range(3):
        a = i*7      #offset =a
        data = [data_set[0 + a],data_set[1 + a],data_set[2 + a],data_set[3 + a],data_set[4 + a],data_set[5 + a],data_set[6 + a]]    
        if i==0:
            b = 0.4
        if i==1:
            b = 0
        if i==2:
            b = -0.4
        bpl = plt.boxplot(data, positions=np.array(range(len(data)))*2.0-b, sym='', widths=0.3)
        if i==0:
            set_box_color(bpl, '#DE2D26') # colors are from http://colorbrewer2.org/
        if i==1:
            set_box_color(bpl, '#0c2c84') # colors are from http://colorbrewer2.org/
        if i==2:
            set_box_color(bpl, '#005824') # colors are from http://colorbrewer2.org/
        plt.xlabel('Data Points per Process')
        plt.ylabel('Time (in seconds)')
        plt.ylim(0.002,100)
        plt.xticks(range(0, len(ticks) * 2, 2), ticks)

    plt.plot([], c='#DE2D26', label='Method1')
    plt.plot([], c='#0c2c84', label='Method2')
    plt.plot([], c='#005824', label='Method3')
    plt.legend()
    if k==0:
        plt.title('Number of Processes = 16')
        plt.savefig('plot16.jpg')
    if k==1:
        plt.title('Number of Processes = 36')
        plt.savefig('plot36.jpg')
    if k==2:
        plt.title('Number of Processes = 49')
        plt.savefig('plot49.jpg')
    if k==3:
        plt.title('Number of Processes = 64')
        plt.savefig('plot64.jpg')   

    plt.show()
