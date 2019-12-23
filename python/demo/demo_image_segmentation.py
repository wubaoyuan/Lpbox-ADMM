import numpy as np
from PIL import Image
import sys
sys.path.append("..")
from functions.util import getUnaryCost, generate_binary_cost, convert2Abc
from functions.lpbox_admm import ADMM_bqp_unconstrained
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix


im = Image.open('cameraman.png')
im = im.convert('L')
node = 5e3
im = im.resize((round(node**0.5), round(node**0.5)), Image.BICUBIC)
I = np.array(im, dtype ='float32') / 255.0

uParams = {'sig':0.1, 'mu_b':0.6, 'mu_f1':0.2, 'mu_f2':0.2}  
unaryCosts = getUnaryCost(I.T.ravel(), uParams)
unaryCosts = np.around(unaryCosts)
W = generate_binary_cost(I)
binary_cost_mat = np.around(3*W)
A,b,c = convert2Abc(unaryCosts,binary_cost_mat)
A = csc_matrix(A)
x0 = np.random.uniform(size = (I.size,1))
x0 = np.array(x0>=0.5, dtype='float32')


ADMM_params = {'std_threshold':1e-6, 'gamma_val':1.0, 'gamma_factor':0.99, \
          'initial_rho':5, 'learning_fact':1+3/100, 'rho_upper_limit':1000, 'history_size':5, 'rho_change_step':5, \
          'rel_tol':1e-5, 'stop_threshold':1e-3, 'max_iters':1e4, 'projection_lp': 2, 'x0':x0}


#optimize
best_sol,x_sol,y1,y2,time_elapsed = ADMM_bqp_unconstrained(A/2, b, ADMM_params)


#show initial image
plt.figure()
image = x0.reshape(I.shape, order='F')
plt.imshow(image, cmap='gray')
plt.imsave('random_initial.png', image, format="png", cmap='gray')


#show result
plt.figure()
image = x_sol.reshape(I.shape, order='F')
plt.imshow(image, cmap='gray')
plt.imsave('segmentation_result.png', image, format="png", cmap='gray')
