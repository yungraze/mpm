# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 11:26:09 2020

@author: Артур
"""

import numpy as np

# Расчет номера элемента прямоугольной эйлеровой сетки,
# которому принадлежит частица

x_min = 1
y_min = 2
nel_x = 5
nel_y = 4
delta_x = 10.0
delta_y = 8.0

xp = 34.5
yp = 25.6

num_x = ((xp - x_min) / delta_x) // 1
num_y = ((yp - y_min) / delta_y) // 1
# 24.5//1 - целая часть числа
print('num_x=', num_x)
print('num_y=', num_y)

x1 = x_min + num_x * delta_x  # Лев.ниж.точка элемента
y1 = y_min + num_y * delta_y

x2 = x1 + delta_x  # Прав.верх.точка лемента с частицей
y2 = y1 + delta_y
print('x1 = ', x1, 'y1=', y1, 'x2=', x2, 'y2=', y2)

N1x = 1 - np.absolute(xp - x1) / delta_x
N2x = 1 - np.absolute(xp - x2) / delta_x
N1y = 1 - np.absolute(yp - y1) / delta_y
N2y = 1 - np.absolute(yp - y2) / delta_y

print('N1x = ', N1x, 'N2x = ', N2x, 'N1y = ', N1y, 'N2y = ', N2y)

N11 = N1x * N1y
N12 = N1x * N2y
N21 = N2x * N1y
N22 = N2x * N2y

print('N11 = ', N11, 'N12 = ', N12, 'N21 = ', N21, 'N22 = ', N22)

# Треугольный элемент сетки (т.е. треугольник АВС)
# и барицентрические координаты в качестве N1, N2, N3

xa = 3  # координаты точки А
ya = 1

xb = 1  # координаты точки В
yb = 2

xc = 1  # координаты точки С
yc = 0

# координаты интересующей нас точки Р
xp = 5 / 3
yp = 1

# Матрица коэффициентов
a11 = xa - xc
a12 = xb - xc
a21 = ya - yc
a22 = yb - yc

# Находим обратную к ней
det = a11 * a22 - a12 * a21
b11 = a22 / det
b12 = -a12 / det
b21 = -a21 / det
b22 = a11 / det

# Находим барицентрич.координаты точки P
N1 = b11 * (xp - xc) + b12 * (yp - yc)
N2 = b21 * (xp - xc) + b22 * (yp - yc)
N3 = 1 - N1 - N2

print('N1=', N1, 'N2=', N2, 'N3=', N3)
