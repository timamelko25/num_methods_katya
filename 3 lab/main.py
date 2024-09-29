import matplotlib.pyplot as plt
import math
import numpy as np


def read_file(file):
    tmp = []
    with open(file, "r") as f:
        tmp = f.read().split()
    return tmp

def read_file2(file):
    tmp = []
    tmp2 = []
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
            tmp2.append(float(data[1].replace(",",".")))
    return tmp, tmp2
    
def research():
    return 0

def main():
    return 0


    
if __name__ == "__main__":
    main()