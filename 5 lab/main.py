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
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
            tmp2.append(float(data[1].replace(",",".")))
            tmp3.append(float(data[2].replace(",",".")))
            tmp4.append(float(data[3].replace(",",".")))
    return tmp, tmp2, tmp3, tmp4
    
def research():
    x, y, dy, err = read_file2("h1_x_y_err.txt")
    x1, y1, dy1, err1 = read_file2("h2_x_y_err.txt")
    
    eps = [pow(10, -i) for i in range(1, 8)]
    
    plt.plot(x1, y1, label="f(x)")
    plt.plot(x, dy, label="h1 = 0.1")
    plt.plot(x1, dy1, label="h2 = 0.01")
    plt.legend()
    plt.title('График функции для 2-ух разибиений')
    plt.xlabel('x') 
    plt.ylabel('y')
    plt.show()
    
    plt.semilogy(x, err, label="err h1 = 0.1")
    plt.semilogy(x1, err1, label="err h2 = 0.01")
    plt.legend()
    plt.title('График ошибки для 2-ух разбиений')
    plt.xlabel('x')
    plt.ylabel('err')
    plt.show()
    
    space = np.linspace(-1, 0, 7)
    
    x, h = read_file("h_eps.txt")
    plt.semilogy(x, h, label="h")
    plt.legend()
    plt.title('График изменения шага по отрезку, eps = 0.0001')
    plt.xlabel('x')
    plt.ylabel('h')
    plt.show()
    
    h, err_eps = read_file("segment-error.txt")
    plt.loglog(eps, err_eps, label="err")
    plt.loglog(eps, eps, label="биссектриса")
    plt.title('График факт ошибки от заданной точности')
    plt.legend()
    plt.xlabel('ошибка')
    plt.ylabel('точность')
    plt.show()
    
    
    
    return 0

def main():
    research()
    return 0


    
if __name__ == "__main__":
    main()