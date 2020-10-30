# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:23:27 2020

@author: Артур
"""

# Выделить все подформулы данной ф.а.в. либо сообщить, что это невозможно

import numpy as np

F = ['(','X','*','(','Y','*','Z',')',')']
#F = ['(','Y','*','Z',')','*','X']
#F = ['X','*','(','Y','*','Z']
#F = ['X']

n = len(F) # Длина (число символов) массива
print('n = ',n)
#print('F = ',F) 

def isLetter(sym): # проверяет. является ли символ одной из букв X, Y, Z
    return (sym=='X')|(sym=='Y')|(sym=='Z')

def isWord(a,b,c,d,e): # проверяет, имеет ли посл-сть вид "( А * В )"
    answer = (a=='(') & (isLetter(b)) & (c=='*') & (isLetter(d)) & (e==')')
    return answer

M = []
count = 0 # счетчик для массива А
A = ['A','B','C','D','E'] # обозначения подформул

isFAV = True
find = False
position = 0
G = [] # заготовка для новой формулы
subs = 'A' # обозначение для заменяемой подформулы

if n==1: # если F состоит из одного символа
    if isLetter(F[0]):
        M.append(F[0])
    else:
        isFAV = False
else:
    if (n < 5):
        isFAV = False    
    else:
        # Найдем вхождение вида "( А * В )"
        i = 0
        while (i < n-5) & (find==False):
            if isWord(F[i],F[i+1],F[i+2],F[i+3],F[i+4]):
                M.append([subs,'=',F[i],F[i+1],F[i+2],F[i+3],F[i+4]])
                find = True
                position = i
            # ===== if ============
            i = i + 1
        # ===== while =======
        if find==True:
            G = []
            for i in range(0,position):
                G.append(F[i])    
            G.append(subs)
            for i in range(position+5,n):
                G.append(F[i])
        if find==False:
            isFAV = False
        # ==== if (find) =======
# === if (len)... else ==============

print('F = ',F)
print('n = ',n)
print('find = ',find)
print('G = ',G)
print('lenG = ',len(G))
    
print('M = ',M)

    

