import calc as c

q_0 = c.stiffness_matrix(0)
q_45 = c.stiffness_matrix(45)
A = c.abd(q_0, q_45)
print(A)

