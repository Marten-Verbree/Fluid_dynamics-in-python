import numpy as np
n = np.zeros((3,3,3))
print(n.size)
def function(z):
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            for k in range(z.shape[2]):
                z[i,j,k] +=1
function(n)
print(n)
