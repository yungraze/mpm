# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:23:27 2020

@author: Артур
"""

# Выделить две подформулы либо сообщить, что это невозможно

import numpy as np

#F = ['X','*','(','Y','*','Z',')']
#F = ['(','Y','*','Z',')','*','X']
F = ['X','*','(','Y','*','Z']
#F = ['X']

n = len(F) # Длина (число символов) массива
print('n = ',n)
print('F = ',F) 

F1 = []
F2 = []

if n == 1:
    print('F -- переменная (не разложима на две формулы)')
else:
    isFAV = True
    # просматриваем первую подформулу
    if (F[0]=='X')|(F[0]=='Y')|(F[0]=='Z'):
        F1.append(F[0])
        i = 1
    else:
        if F[0]=='(':
            count = 0
            i = 0
            j = i
            while (j == i)|((count > 0)&(j < n)):
                if F[j]=='(':
                    count = count + 1
                if F[j]==')':
                    count = count - 1
                F1.append(F[j])
                j = j + 1
            #===== while ===========
            i = j
            if (count < 0)|(count > 0):
                isFAV = False
        else:
            isFAV = False
            print(F[0],'-- ни переменная ни открыв. скобка')
        #====== if else =================
    #======= if else ====================            

    if (F[i]=='*')==False:
        isFAV = False
        print('После F1 стоит символ, отличный от *')
    i = i + 1
    
    # просматриваем вторую подформулу
    if (F[i]=='X')|(F[i]=='Y')|(F[i]=='Z'):
        F2.append(F[i])
    else:
        if F[i]=='(':
            count = 0
            j = i
            while (j == i)|((count > 0)&(j < n)):
                if F[j]=='(':
                    count = count + 1
                if F[j]==')':
                    count = count - 1
                F2.append(F[j])
                j = j + 1
            #===== while ===========
            i = j
            if (count < 0)|(count > 0):
                isFAV = False
        else:
            isFAV = False
        # ====if F[i]=(... else =========
    # ===== if (F[i]==X)... else ==========
    print('i = ',i) 
    if i < n - 1:
        isFAV = False

# ==== if else ============
        
if (isFAV==True):
    print('F1 = ',F1)
    print('F2 = ',F2)
else:
    print('F -- не ф.а.в.')

    
    

    

