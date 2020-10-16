# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:23:27 2020

@author: Артур
"""

# Проверка: образует ли посл-сть символов формулу алгебры высказываний

import numpy as np

F = ['(',')','(',')',')','(','(',')']

n = len(F) # Длина (число символов) массива

print(F) 
   
isFAV = True
count = 0

# Просматриваем весь массив F
for i in range(0,n):
    if F[i]=='(':
        count = count + 1
    if F[i]==')':
        count = count - 1
    if count < 0 : # в какой-то момент (число закр.скобок) > (числа откр.)
        isFAV = False 

# 1) Ни в какой момент число закр. не должно быть > откр (isFAV = True))
# 2) Общее число откр должно быть = числу закр. (count = 0)
if (isFAV==True)&(count==0):
    print('Это ф.а.в.')
else:
    print('Это не ф.а.в.')
    
# 
    

    

