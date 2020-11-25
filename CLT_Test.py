import numpy as np
import matplotlib.pyplot as plt

def CLT_Shaft(stack):
    minSafety = np.sum(stack[0]*stack[1])

    return 1/minSafety