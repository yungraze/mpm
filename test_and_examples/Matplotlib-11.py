# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 11:54:12 2020

@author: Артур
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

fig = plt.figure()

plt.scatter(0.0, 0.0)
plt.plot([0, 1, 0, 2], [1, 0, 1, 3])

plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.show()
