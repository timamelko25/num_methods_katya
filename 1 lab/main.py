import matplotlib.pyplot as plt
import math
import numpy as np

def func1(x):
    return (2*x-math.log10(x))

def sgn(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0
    
def func2(x):
    return sgn(x)*pow(x,4)-18*pow(x,2) + 6

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
    
def research(filename, func, a, b):
    x_5, y_5 = read_file2(f"{filename}_mesh_5.txt")
    x_7, y_7 = read_file2(f"{filename}_mesh_7.txt")
    x_10, y_10 = read_file2(f"{filename}_mesh_10.txt")
    
    x_5_10000, y_5_10000 = read_file2(f"{filename}_mesh_5_10000.txt")
    x_7_10000, y_7_10000 = read_file2(f"{filename}_mesh_7_10000.txt")
    x_10_10000, y_10_10000 = read_file2(f"{filename}_mesh_10_10000.txt")
    
    x2_5, y2_5 = read_file2(f"{filename}_chebish_5.txt")
    x2_7, y2_7 = read_file2(f"{filename}_chebish_7.txt")
    x2_10, y2_10 = read_file2(f"{filename}_chebish_10.txt")
    
    x2_5_10000, y2_5_10000 = read_file2(f"{filename}_chebish_5_10000.txt")
    x2_7_10000, y2_7_10000 = read_file2(f"{filename}_chebish_7_10000.txt")
    x2_10_10000, y2_10_10000 = read_file2(f"{filename}_chebish_10_10000.txt")
    
    space = np.linspace(a, b, 10000)
    f_x = [func(x) for x in space]
    
    plt.plot(x_5, y_5, label="5", marker='o', markersize=7)
    plt.plot(x2_5, y2_5, label="5", marker='o', markersize=7)
    plt.plot(space, f_x, label="5_10000")
    plt.plot(x_5_10000, y_5_10000, label="5_10000")
    plt.plot(x2_5_10000, y2_5_10000, label="5_10000")
    plt.title('Первый график') 
    plt.legend() # дописать
    plt.show()
    
    plt.plot(x_7, y_7, label="7", marker='o', markersize=7)
    plt.plot(x2_7, y2_7, label="7", marker='o', markersize=7)
    plt.plot(space, f_x, label="7_10000")
    plt.plot(x_7_10000, y_7_10000, label="7_10000")
    plt.plot(x2_7_10000, y2_7_10000, label="7_10000")
    plt.title('Первый график') 
    plt.legend() # дописать
    plt.show()
    
    plt.plot(x_10, y_10, label="10", marker='o', markersize=7)
    plt.plot(x2_10, y2_10, label="10", marker='o', markersize=7)
    plt.plot(space, f_x, label="10_10000")
    plt.plot(x_10_10000, f_x, label="10_10000")
    plt.plot(x_10_10000, y_10_10000, label="10_10000")
    plt.plot(x2_10_10000, y2_10_10000, label="10_10000")
    plt.title('Первый график') 
    plt.legend() # дописать
    plt.show()
    
    x, error_5 = read_file2(f"{filename}_mesh_error_5.txt")
    x, error_7 = read_file2(f"{filename}_mesh_error_7.txt")
    x, error_10 = read_file2(f"{filename}_mesh_error_10.txt")
    
    plt.semilogy(x, error_5, label="5")
    plt.title('Первый график') 
    plt.legend() # дописать
    plt.show()
    
    plt.semilogy(x, error_7, label="7")
    plt.title('Первый график') 
    plt.legend() # дописать
    plt.show()
    
    plt.semilogy(x, error_10, label="10")
    plt.title('Первый график') 
    plt.legend() # дописать
    plt.show()
    
    max_error_mesh = read_file(f"{filename}_mesh_max_error.txt")
    max_error_chebish = read_file(f"{filename}_chebish_max_error.txt")
    n = [i+5 for i in range(96)]
    plt.semilogy(n, max_error_mesh, label="5")
    plt.semilogy(n, max_error_chebish, label="5")
    plt.title('Первый график')
    plt.legend() # дописать
    plt.show()
    
    point1_mesh = read_file(f"{filename}_mesh_point1.txt")
    point1_chebish = read_file(f"{filename}_chebish_point1.txt")
    plt.semilogy(n, point1_mesh, label="5")
    plt.semilogy(n, point1_chebish, label="5")
    plt.title('Первый график')
    plt.legend() # дописать
    plt.show()
    
    point2_mesh = read_file(f"{filename}_mesh_point2.txt")
    point2_chebish = read_file(f"{filename}_chebish_point2.txt")
    plt.semilogy(n, point2_mesh, label="5")
    plt.semilogy(n, point2_chebish, label="5")
    plt.title('Первый график')
    plt.legend() # дописать
    plt.show()
    
    

def main():
    # добавить подписи осей
    research("func1", func1, 1, 10)
    research("func2", func2, 0, 2)
    return 0


    
if __name__ == "__main__":
    main()