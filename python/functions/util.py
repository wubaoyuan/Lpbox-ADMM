import numpy as np
import math


def getUnaryCost(nodes, params): 
    sig = params['sig']
    mu_b = params['mu_b']
    mu_f1 = params['mu_f1']
    mu_f2 = params['mu_f2']    

    const = 1./2*np.log(2*math.pi) + np.log(sig)

    alpha_b = (nodes - mu_b)**2/(2*(sig**2)) + const
    
    aa = np.exp(-(nodes - mu_f1)**2/(2*(sig**2))) +  np.exp(-(nodes - mu_f2)**2/(2*(sig**2)))
    alpha_f = -np.log(aa + 2.2204e-16) + const + np.log(2)
    
    cost = np.vstack((alpha_b,alpha_f))

    return cost



def generate_binary_cost(I):
    sigma = np.std(I)
    neighborhood_edge = 1

    nrows,ncols = I.shape
    pairs = generate_pixel_pairs(nrows, ncols, neighborhood_edge)

    I_new = I.T.ravel()
    
    diff_vec = ((I_new[pairs[:,0]]-I_new[pairs[:,1]])**2)/sigma
    diff_vec = np.exp(-diff_vec)
    
    binary_cost_mat = np.zeros((I.size, I.size))
    binary_cost_mat[pairs[:,0],pairs[:,1]] = diff_vec
    
    return binary_cost_mat




def generate_pixel_pairs(nrows,ncols,k):

    offsets = get_offsets(k)
    
    rmat,cmat = np.meshgrid(np.array(range(1,nrows+1)), np.array(range(1,ncols+1)))
    spatial_indices = np.hstack((rmat.reshape(-1,1),cmat.reshape(-1,1)))
    
    num_pixels = rmat.size

    for p in range(offsets.shape[0]):
        offset = offsets[p,:];
        if np.linalg.norm(offset)>0:
            temp1 = np.tile(offset, (num_pixels, 1))
            temp2 = np.hstack((spatial_indices, spatial_indices+temp1))
            
            if p==0:
                neighborhood_mat = temp2
            else:
                neighborhood_mat = np.vstack((neighborhood_mat, temp2))


    valid_idx = (neighborhood_mat[:,2]>=1) & (neighborhood_mat[:,2]<=nrows) & (neighborhood_mat[:,3]>=1) & (neighborhood_mat[:,3]<=ncols)
    
    coor1 = np.vstack((neighborhood_mat[valid_idx,0].reshape(1,-1), neighborhood_mat[valid_idx,1].reshape(1,-1)))
    ind1 = np.ravel_multi_index(coor1-1, (nrows,ncols), order='F')
    
    coor2 = np.vstack((neighborhood_mat[valid_idx,2].reshape(1,-1), neighborhood_mat[valid_idx,3].reshape(1,-1)))
    ind2 = np.ravel_multi_index(coor2-1, (nrows,ncols), order='F')
    pairs = np.hstack((ind1.reshape(-1,1),ind2.reshape(-1,1)))


    return pairs
    
    
    

def get_offsets(k):
    col_offsets = np.tile(np.array(range(-k,k+1)), (2*k+1,1))
    row_offsets = col_offsets.T
    offsets = np.vstack((col_offsets.reshape(-1),row_offsets.reshape(-1))).T

    return offsets 
    

def convert2Abc(unary, pairwise):
    b = unary[1,:]-unary[0,:]
    b=b.reshape(-1,1)
    
    c = sum(unary[0,:])
    
    A  = -pairwise
    We = -np.sum(A,1)
    
    A = A+ np.diag(We)
    A = 2*A
    
    return A,b,c
