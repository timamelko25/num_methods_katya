import matplotlib.pyplot as plt
import math
import numpy as np

def f(x):
    return math.sin(x)

def read_file(file):
    tmp = []
    tmp2 = []
    tmp3 = []
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
            tmp2.append(float(data[1].replace(",",".")))
            tmp3.append(float(data[1].replace(",",".")))
    return tmp, tmp2, tmp3

def read_file2(file):
    tmp = []
    with open(file, "r") as f:
        for line in f:
            data = line.split()
            tmp.append(float(data[0].replace(",",".")))
    return tmp
    
def research():
    x1, y1, err1 = read_file("h1_x_y_err.txt")
    x2, y2, err2 = read_file("h2_x_y_err.txt")
    x = np.linspace(0, math.pi/2, 100)
    f_x = [f(x) for x in x]
    
    plt.plot(x1, y1, label="h1 = 0.16")
    plt.plot(x2, y2, label="h2 = 0.016")
    plt.plot(x, f_x, label="f(x)")
    plt.legend()
    plt.title('График функции для 2-ух разибиений')
    plt.xlabel('x') 
    plt.ylabel('y')
    plt.show()
    
    
    eps = [pow(10, -i) for i in range(1, 7)]
    
    plt.semilogy(x1, err1, label="err h1 = 0.16")
    plt.semilogy(x2, err2, label="err h2 = 0.016")
    plt.legend()
    plt.title('График ошибки для 2-ух разбиений')
    plt.xlabel('x')
    plt.ylabel('err')
    plt.show()
    
    
    err3 = read_file2("error_eps.txt")
    h = read_file2("runge_h_2.txt")
    
    x = np.linspace(0, math.pi/2, 1034)
    
    plt.semilogy(x, h, label="h")
    plt.legend()
    plt.title('График изменения шага по отрезку, eps = 0.01')
    plt.xlabel('x')
    plt.ylabel('h')
    plt.show()
    
    plt.loglog(eps, err3, label="err")
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