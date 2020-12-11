import numpy as np
from PIL import Image
import sys, getopt
sys.path.append("../python/")
from functions.util import getUnaryCost, generate_binary_cost, convert2Abc
from functions.lpbox_admm import ADMM_bqp_unconstrained
import matplotlib.pyplot as plt
import itertools
from scipy.sparse import csc_matrix
from _CInterface import ffi, lib

argv = sys.argv
recvfd = int(argv[1])
sendfd = int(argv[2])
inputfile = argv[3]
outputfile = argv[4]

im = Image.open(inputfile).convert('L')

im.save('out/orig_' + inputfile, format="png")
# print(im.size)
node = 5e3
im = im.resize((round(node**0.5), round(node**0.5)), Image.BICUBIC)
I = np.array(im, dtype ='float64') / 255.0

uParams = {'sig':0.1, 'mu_b':0.6, 'mu_f1':0.2, 'mu_f2':0.2}
unaryCosts = getUnaryCost(I.T.ravel(), uParams)
unaryCosts = np.around(unaryCosts)
W = generate_binary_cost(I)
binary_cost_mat = np.around(3*W)
A,b,c = convert2Abc(unaryCosts,binary_cost_mat)
A = csc_matrix(A)
# x0 = loadtxt('x0', dtype='float64')
# x0 = x0.reshape((x0.shape[0], 1))
# print(x0.shape)
x0 = np.random.uniform(size = (I.size,1))
x0 = np.array(x0>=0.5, dtype='float64')
# print(x0.shape)

vecSize = len(b)
print("b size: ", vecSize)
print("A size: ", A.getnnz())

rfd = ffi.cast("int", recvfd)
bLen = ffi.cast("int", vecSize)
_bLen = ffi.new("int[1]", [vecSize])
_b = ffi.cast("void *", b.ctypes.data)
_x0 = ffi.cast("void *", x0.ctypes.data)
sfd = ffi.cast("int", sendfd)
ALen = ffi.cast("int", A.getnnz())
_ALen = ffi.new("int[1]", [A.getnnz()])

lib.write(sfd, _bLen, ffi.cast("int", 4))
lib.write(sfd, _x0, ffi.cast("int", 8 * vecSize))
lib.write(sfd, _b, ffi.cast("int", 8 * vecSize))
lib.write(sfd, _ALen, ffi.cast("int", 4))

Ax = A.tocoo()
for i, j, v, in zip(Ax.row, Ax.col, Ax.data):
    # print(i, j, v)
    _i = ffi.cast("int", i)
    _j = ffi.cast("int", j)
    _v = ffi.cast("double", v)
    lib.sendTriplet(_i, _j, _v, sfd)

_x_sol = ffi.new("double[]", vecSize)
lib.read(rfd, _x_sol, ffi.cast("int", 8 * vecSize))
print("read x_sol from c algorithm")

tmp = ffi.buffer(_x_sol, 8 * vecSize)
x_sol = np.frombuffer(tmp, dtype=np.float_, count=vecSize)

print(x_sol)

plt.figure()
image = x0.reshape(I.shape[0:2], order='F')
plt.imshow(image, cmap='gray')
plt.imsave('out/init_' + inputfile, image, format="png", cmap='gray')

#show result
plt.figure()
image = x_sol.reshape(I.shape[0:2], order='F')
plt.imshow(image, cmap='gray')
plt.imsave('out/seg_' + inputfile, image, format="png", cmap='gray')

print("Excecute success")