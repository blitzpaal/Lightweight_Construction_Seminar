import numpy as np

d = np.ones(3)
a = np.reshape(d,(3,1))
b = np.identity(3)
print(a)
print(b)

c = np.concatenate((b,a),axis=1)

print(c)