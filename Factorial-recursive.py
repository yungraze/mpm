# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:23:27 2020

@author: Артур
"""

# Вычисление факториала с помощью рекурсии

import numpy as np

def Factorial(n):
    if n > 1:
        return n * Factorial(n-1)
    else:
        return 1

n = 3
fact = Factorial(n)
print('fact = ',fact)
