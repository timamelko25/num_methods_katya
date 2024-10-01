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
    tmp3 = []
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
            tmp2.append(float(data[1].replace(",",".")))
            tmp3.append(float(data[2].replace(",",".")))
    return tmp, tmp2, tmp3
    
def research(filename):
    i, err, h = read_file2(filename)
    eps = [pow(10, -i) for i in range(1, 8)]
    plt.loglog(err, eps, label="err")
    plt.loglog(eps, eps, label="биссектриса")
    plt.title('График факт ошибки от заданной точности')
    plt.legend()
    plt.xlabel('ошибка')
    plt.ylabel('точность')
    plt.show()
    
    plt.semilogx(eps, i, label="iter")
    plt.title('График итераций от заданной точности')
    plt.legend()
    plt.xlabel('точность')
    plt.ylabel('итерации')
    plt.show()
    
    for j in range(7):
        h[j] = math.log2(h[j])
    
    plt.semilogx(err, h, label="log2(h)")
    plt.title('График факт ошибки от длины разбиения')
    plt.legend()
    plt.xlabel('ошибка')
    plt.ylabel('длина разбиения')
    plt.show()
    return 0

def main():
    research("data.txt")
    return 0


    
if __name__ == "__main__":
    main()