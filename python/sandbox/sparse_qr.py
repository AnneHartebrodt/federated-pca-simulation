
import numpy as np
import numpy.linalg as na
import scipy.linalg as la
mat = [[1.2,3.2,0,0,0, 0, 0, 1.2],[0,0,2.5, 0,0, 0, 0, 0.2],[0,0,0,0,0, 1.1, 0, 1.9],[0,0,3.3,0,1.15, 0, 0, 0]]
mat = np.matrix(mat.T)

b = [3,2,5,1,2,3,1,0,2]

q,r= la.qr(mat)