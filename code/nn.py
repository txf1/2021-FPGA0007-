import numpy
import scipy.io as scipy
import math

train_label_normal = '/home/xilinx/jupyter_notebooks/DDC/train_label_normal.mat'
test_label_normal = '/home/xilinx/jupyter_notebooks/DDC/test_label_normal.mat'
test_normal = '/home/xilinx/jupyter_notebooks/DDC/test_normal.mat'
train_normal = '/home/xilinx/jupyter_notebooks/DDC/train.mat'

test = '/home/xilinx/jupyter_notebooks/DDC/test.mat'
test_label = '/home/xilinx/jupyter_notebooks/DDC/test_label.mat'
train = '/home/xilinx/jupyter_notebooks/DDC/train.mat'
train_label = '/home/xilinx/jupyter_notebooks/DDC/train_label.mat'

train_label_normal = scipy.loadmat(train_label_normal)
test_label_normal = scipy.loadmat(test_label_normal)
test_normal = scipy.loadmat(test_normal)
train_normal = scipy.loadmat(train_normal)
test = scipy.loadmat(test)
test_label = scipy.loadmat(test_label)
train = scipy.loadmat(train)
train_label = scipy.loadmat(train_label)

train_label_normal = train_label_normal['train_label_normal']
test_label_normal = test_label_normal['test_label_normal']
test_normal = test_normal['test_normal']
train_normal = train_normal['train']
test = test['test']
test_label = test_label['test_label']
train = train['train']
train_label = train_label['train_label']

for itr in train_normal:
	itr = itr[0:292]

for itr in test_normal:
	itr = itr[0:292]

train = list(train)
train_label = list(train_label)
train_normal = list(train_normal)
train_label_normal = list(train_label_normal)

test = list(test)
test_label = list(test_label)
test_normal = list(test_normal)
test_label_normal = list(test_label_normal)

train.extend(train_normal)
train_label.extend(train_label_normal)

test.extend(test_normal)
test_label.extend(test_label_normal)

syn0 = 2*numpy.random.random((293,10)) - 1
syn1 = 2*numpy.random.random((10,1)) - 1
for j in range(60000):
    l1 = 1/(1+numpy.exp(-(numpy.dot(train,syn0))))
    l2 = 1/(1+numpy.exp(-(numpy.dot(l1,syn1))))
    l2_delta = (train_label - l2)*(l2*(1-l2))
    l1_delta = l2_delta.dot(syn1.T) * (l1 * (1-l1))
    syn1 += l1.T.dot(l2_delta)
    syn0 += X.T.dot(l1_delta)