import numpy as np

A = np.ones((3,3))
B = np.zeros((3,3))

D = np.ones((3,3))*2



# c = np.concatenate((b,a),axis=1)


ABD = np.vstack((np.hstack((A,B)),np.hstack((B,D))))
e = np.ones(5)
print(e)
f = e[0:3]
print(f)
