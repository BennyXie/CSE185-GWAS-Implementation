# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import csv
import os
import sys
import numpy as np
import scipy.stats as stats
def read_phen(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        phen = []
        for row in reader:
            phen.append([str(row[1]),float(row[2])])
    return phen


def read_gen(filename):
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        row = next(reader)
        while row[0] != "#CHROM":
            row = next(reader)
        keys = row
        sampleStart=keys.index("FORMAT")+1

        gen = {}
        for row in reader:
            if row[0][0] != "#":
                gen[str(row[2])] = {}
                for i in range(sampleStart, len(row)):
                    genotype=-1
                    if row[i] == "0/0":
                        genotype = 0
                    elif row[i] == "0/1":
                        genotype = 1
                    elif row[i] == "1/1":
                        genotype = 2
                    gen[str(row[2])][str(keys[i])] = genotype
    return gen

def regression(x, y):
    slope, intercept, r, p, se = stats.linregress(x, y)
    return slope, p
# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    gen = read_gen("/var/home/yifeiding/cse185/lab3/gtdata_1000Genomes_pruned.vcf")
    phen = read_phen("/var/home/yifeiding/public/lab3/lab3_gwas.phen")

    for i in gen.keys():
        x = []
        y = []
        for j in range(len(phen)):
            #if gen[i][phen[j][0]] != -1:
            x.append(gen[i][phen[j][0]])
            y.append(phen[j][1])
        slope, p = regression(x, y) # This part doesn't work
        print(i, p)
