#!/usr/bin/python3.8
#coding: utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def parse(line):
    if "buffer" not in line:
        return -1, -1
    strs = line.split(':')
    return int(strs[0].split(' ')[-1].strip()), int(strs[-1].strip())
        


if __name__ == "__main__":
    list = []
    
    with open(sys.argv[1]) as f:
        while True:
            line = f.readline().strip('\n')
            if line == "":
                break
            x, y  = parse(line)
            if x == -1 or y == -1:
                continue
            if x >= len(list):
                list.append([])
            else:
                list[x].append(y)            
    
    xaxis = np.arange(0, len(list[0]), 1)
    
    
    for thread in list:
        plt.plot(xaxis, thread)
    plt.savefig('./plot.png')