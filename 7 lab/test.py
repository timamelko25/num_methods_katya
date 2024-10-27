import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Параметры задачи
a, b = 0, np.pi
A = 0
B = 1

# Уравнение: y'' + y = 2x - pi
# Сначала перепишем его как систему двух ОДУ первого порядка:
# z1 = y, z2 = y'
# z1' = z2
# z2' = 2x - pi - z1

def ode_system(x, z):
    z1, z2 = z
    dz1dx = z2
    dz2dx = 2 * x - np.pi - z1
    return [dz1dx, dz2dx]

# Решение первой вспомогательной задачи для y1, где y1(0) = 0, y1'(0) = 1
sol1 = solve_ivp(ode_system, [a, b], [0, 1], t_eval=np.linspace(a, b, 100))

# Решение второй вспомогательной задачи для y2, где y2(0) = 1, y2'(0) = 0
sol2 = solve_ivp(ode_system, [a, b], [1, 0], t_eval=np.linspace(a, b, 100))

# Определение коэффициентов суперпозиции C1 и C2
y1_b, y1p_b = sol1.y[0, -1], sol1.y[1, -1]
y2_b, y2p_b = sol2.y[0, -1], sol2.y[1, -1]

C1 = (B - y1_b) / (y2_b + y2p_b)
C2 = -C1 * y1p_b

# Финальное решение как линейная комбинация y = C1 * y1 + C2 * y2
y = C1 * sol1.y[0] + C2 * sol2.y[0]

# Построение графика численного решения
x_vals = np.linspace(a, b, 100)
exact_solution = 2 * x_vals - np.pi + np.pi * np.cos(x_vals) + np.sin(x_vals)

plt.plot(x_vals, y, label='Numerical Solution')
plt.plot(x_vals, exact_solution, '--', label='Exact Solution')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.title("Numerical Solution vs Exact Solution")
plt.grid()
plt.show()
