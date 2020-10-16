# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:23:27 2020

@author: Артур
"""

# Выделить все подформулы данной ф.а.в. либо сообщить, что это невозможно

import numpy as np

F = ['(','(','X','*','Y',')','*','(','X','*','Z',')',')']
#F = ['(','X','*','(','Y','*','Z',')',')']
#F = ['(','(','Y','*','Z',')','*','X',')']
#F = ['X','*','(','Y','*','Z']
#F = ['X']

n = len(F) # Длина (число символов) массива
print('n = ',n)
#print('F = ',F) 

def isLetter(sym): # проверяет. является ли символ одной из букв X, Y, Z
    return (sym=='X')|(sym=='Y')|(sym=='Z')|(sym=='A')|(sym=='B')

def isWord(a,b,c,d,e): # проверяет, имеет ли посл-сть вид "( А * В )"
    answer = (a=='(') & (isLetter(b)) & (c=='*') & (isLetter(d)) & (e==')')
    return answer

Subformulas = []
count = 0 # счетчик для массива А
A = ['A','B','C','D','E'] # обозначения подформул

isFAV = True
#find = False
position = 0
G = [] # заготовка для новой формулы
subs = 'A' # обозначение для заменяемой подформулы

def Simplify(F,subs):
    isFAV = True
    find = False
    position = 0
    variable = False
    G = []
    M = []
    n = len(F)
    if n==1: # если F состоит из одного символа
        if isLetter(F[0]):
            G.append(F[0])
            variable = True
        else:
            isFAV = False
    else:
        if (n < 5):
            isFAV = False    
        else:
            # Найдем вхождение вида "( А * В )"
            i = 0
            while (i < n-4) & (find==False):
                if isWord(F[i],F[i+1],F[i+2],F[i+3],F[i+4]):
                    M = [subs,'=',F[i],F[i+1],F[i+2],F[i+3],F[i+4]]
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
        # == if (n<5)... else ==============
    # === if n==1... else ==============
    return [G, M, isFAV,variable]
# ===== def Simplify =============

print('Исходная формула:')
print('F = ',F)
length = len(F)
isFAV_F = True
if (length==1)&(isLetter(F[0])==False):
    isFAV_F==False
count = 0
t = 0
while (length > 1)&(t < 10)&(isFAV_F==True) :
    S = Simplify(F,A[count])
    count = count + 1
    print('Шаг',t)
    print('Выделенная подформула:',S[1])
    print('Упрощенная формула:',S[0])
    isFAV_F = S[2]
    length = len(S[0])
    F = []
    for i in range(0,len(S[0])):
        F.append(S[0][i])
    t = t + 1
    
# ====== while ===========

print()
if isFAV_F == True:
    print('F -- формула')    
if isFAV_F==False:
    print('F -- не формула')

