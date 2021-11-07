#show_pkl.py

import pickle

path = '/home/xilinx/jupyter_notebooks/IACM-master-main/output/mit-bih-database.pkl'

f = open(path,'rb')
data = pickle.load(f)

print(data)
