import numpy as np
import scipy as sp
from scipy.sparse import csc_matrix, coo_matrix, csr_matrix, lil_matrix, random
from jinja2 import Environment
from jinja2.loaders import FileSystemLoader
import os

seed = 42
density = 0.1
gamma = 0.01

out_file_name = 'qp_data.hpp'
in_file_name  = 'qp_data.in.hpp'

NV = 100
NC = 10

H = csc_matrix((NV, NV))
A = csc_matrix((NV, NV))

myinf = 1e10

for i in range(NV):
    H[i,i] = 1.0

# H = H + gamma*random(NV, NV, density=density, format='csc', random_state=seed)
# H = H.T*H

for i in range(NC):
    A[i,i] = 1.0

H_ri =  H.indices 
H_cp =  H.indptr 
H_val = H.data 
H_nnz = H.nnz

A_ri =  A.indices 
A_cp =  A.indptr 
A_val = A.data 
A_nnz = A.nnz

g = np.ones((NV,1))

lb = -myinf*np.ones((NV, 1))
ub = myinf*np.ones((NV, 1))

lbA = -np.ones((NC, 1))
ubA = np.ones((NC, 1))

print('rendering templated C++ code...')
env = Environment(loader=FileSystemLoader(os.path.dirname(os.path.abspath(__file__))))
tmpl = env.get_template(in_file_name)

code = tmpl.render(NV = NV, NC = NC, H_cp = H_cp, H_ri = H_ri, H_val = H_val, H_nnz = H_nnz, \
        A_cp = A_cp, A_ri = A_ri, A_val = A_val, A_nnz = A_nnz, g = g, lb = lb, ub = ub, lbA = lbA, ubA = ubA)

with open(out_file_name, "w+") as f:
    f.write(code.replace('inf', 'Inf'))
