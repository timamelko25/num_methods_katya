import matplotlib.pyplot as plt
import math
import numpy as np


def read_file(file):
    tmp = []
    tmp2 = []
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
            tmp2.append(float(data[1].replace(",",".")))
    return tmp, tmp2

def read_file2(file):
    tmp = []
    tmp2 = []
    tmp3 = []
    tmp4 = []
    tmp5 = []
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
            tmp2.append(float(data[1].replace(",",".")))
            tmp3.append(float(data[2].replace(",",".")))
            tmp4.append(float(data[3].replace(",",".")))
            tmp5.append(float(data[4].replace(",",".")))
    return tmp, tmp2, tmp3, tmp4, tmp5
    
def research():
    h, x1, err1, x2, err2 = read_file2("data.txt")
    
    plt.plot(h, err1, label="err")
    plt.legend()
    plt.grid(True)
    plt.title('График ошибки для n разибиений Метод Рунге-Кутта')
    plt.xlabel('h') 
    plt.ylabel('err')
    plt.show()
    
    plt.plot(h, x1, label="h")
    plt.legend()
    plt.grid(True)
    plt.title('График координаты для n разибиений Метод Рунге-Кутта')
    plt.xlabel('h') 
    plt.ylabel('x')
    plt.show()
    
    plt.plot(h, err2, label="err")
    plt.legend()
    plt.grid(True)
    plt.title('График ошибки для n разибиений Метод Предиктор Корректор')
    plt.xlabel('h') 
    plt.ylabel('err')
    plt.show()
    
    plt.plot(h, x2, label="h")
    plt.legend()
    plt.grid(True)
    plt.title('График координаты для n разибиений Метод Предиктор Корректор')
    plt.xlabel('h') 
    plt.ylabel('x')
    plt.show()

def main():
    research()

if __name__ == "__main__":
    main()