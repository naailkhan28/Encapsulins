import matplotlib.pyplot as plt
import numpy as np

center1 = (50, 60)
center2 = (80, 20)
center3 = (90, 120)
distance = 20


x1 = np.random.uniform(center1[0], center1[0] + distance, size=(100,))
y1 = np.random.normal(center1[1], distance, size=(100,)) 

x2 = np.random.uniform(center2[0], center2[0] + distance, size=(100,))
y2 = np.random.normal(center2[1], distance, size=(100,)) 

x3 = np.random.uniform(center3[0], center3[0] + distance, size=(100,))
y3 = np.random.normal(center3[1], distance, size=(100,)) 

fig,ax = plt.subplots(1)

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.scatter(x1, y1)
ax.scatter(x2, y2)
ax.scatter(x3, y3)
plt.show()