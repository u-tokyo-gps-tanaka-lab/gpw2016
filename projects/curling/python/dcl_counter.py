# -*- coding: utf-8 -*-
# dcl_counter.py
# Katsuki Ohto

import argparse
import glob
import numpy as np

def analyze_logs(logs):
    # result survey
    wl = {}
    scores = {}
    # simulation survey
    errors = [[], []]
    
    for f in glob.glob(logs):
        # get data
        name = ["", ""]
        chosen = (0, 0)
        run = (0, 0)
        score = [-1, -1]
        flip = 1
        
        for line in open(f, 'r'):
            data = line.split()
            if 'First=' in data[0]:
                name[0] = data[0][6:]
            elif 'Second=' in data[0]:
                name[1] = data[0][7:]
            elif 'BESTSHOT' in data[0]:
                chosen = (float(data[1]), float(data[2]))
            elif 'RUNSHOT' in data[0]:
                run = (float(data[1]), float(data[2]))
                errors[0].append(run[0] - chosen[0])
                errors[1].append(run[1] - chosen[1])
            elif 'TOTALSCORE' in data[0]:
                score[0] = int(data[1])
                score[1] = int(data[2])
            elif 'SCORE' in data[0]:
                if flip * int(data[1]) < 0:
                    flip = -flip
    
        for i in range(2):
            if name[i] not in wl:
                wl[name[i]] = [[0, 0, 0, 0], [0, 0, 0, 0]]
                scores[name[i]] = [[0, 0], [0, 0]]

        for i in range(2):
            for c in range(2):
                scores[name[i]][i][i ^ c] += score[c]

        if score[0] > score[1]:
            wl[name[0]][0][0] += 1
            wl[name[1]][1][3] += 1
        elif score[0] < score[1]:
            wl[name[0]][0][3] += 1
            wl[name[1]][1][0] += 1
        else:
            if flip == -1:
                wl[name[0]][0][1] += 1
                wl[name[1]][1][2] += 1
            else:
                wl[name[0]][0][2] += 1
                wl[name[1]][1][1] += 1
    
    print(wl)
    print(scores)
    print("error in Vx : mean = %f stddev = %f" % (np.mean(errors[0]), np.std(errors[0])))
    print("error in Vy : mean = %f stddev = %f" % (np.mean(errors[1]), np.std(errors[1])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--logs', required=True)
    args = parser.parse_args()
    
    analyze_logs(args.logs)


