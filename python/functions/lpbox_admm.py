import numpy as np
import time
import scipy.sparse.linalg as linalg
from scipy.sparse import csc_matrix


"""
min_x (x^T A x + b^T x) such that x is {0,1}^n
ADMM update steps with x0 being feasible and binary
"""
def ADMM_bqp_unconstrained(A, b, all_params=None):
    initial_params = {'std_threshold':1e-6, 'gamma_val':1.0, 'gamma_factor':0.99, \
          'initial_rho':5, 'learning_fact':1+3/100, 'rho_upper_limit':1000, 'history_size':5, 'rho_change_step':5, \
          'rel_tol':1e-5, 'stop_threshold':1e-3, 'max_iters':1e4, 'projection_lp': 2, 'pcg_tol':1e-3, 'pcg_maxiters':1e3}
    
    if all_params==None:
        all_params = initial_params
    else:
        for k in initial_params.keys():
            if k not in all_params.keys():
                all_params[k] = initial_params[k]
    
    n = b.size
    stop_threshold = all_params['stop_threshold']
    std_threshold = all_params['std_threshold']
    max_iters = all_params['max_iters']
    initial_rho = all_params['initial_rho']
    rho_change_step = all_params['rho_change_step']
    gamma_val = all_params['gamma_val']
    learning_fact = all_params['learning_fact']
    history_size = all_params['history_size']
    projection_lp = all_params['projection_lp']
    gamma_factor = all_params['gamma_factor']
    pcg_tol = all_params['pcg_tol']
    pcg_maxiters = all_params['pcg_maxiters']
    
    
    x_sol = all_params['x0']
    y1 = x_sol
    y2 = x_sol
    z1 = np.zeros_like(y1)
    z2 = np.zeros_like(y2)
    rho1 = initial_rho
    rho2 = rho1
    obj_list = []
    std_obj = 1
    
    
    # initiate the binary solution
    prev_idx = x_sol
    best_sol = prev_idx
    best_bin_obj = compute_cost(best_sol,A,b)
    
    
    time_elapsed = 0
    for iter in range(int(max_iters)):
        t1 = time.time()
        
        # update y1
        y1 = project_box(x_sol+z1/rho1)
    
        # update y2
        y2 = project_shifted_Lp_ball(x_sol+z2/rho2, projection_lp)
        
        row = np.array(range(n))
        colum = np.array(range(n))
        data = (rho1+rho2)*np.ones(n)
        sparse_matrix = csc_matrix((data, (row, colum)), shape=(n, n))
        x_sol,cg_flag = linalg.cg(2*A+sparse_matrix, rho1*y1+rho2*y2-(b+z1+z2), y1, tol=pcg_tol, maxiter=pcg_maxiters)
        x_sol = x_sol.reshape(-1,1)

        
        # update z1 and z2
        z1 = z1+gamma_val*rho1*(x_sol-y1)
        z2 = z2+gamma_val*rho2*(x_sol-y2)

        
        t2 = time.time()
        time_elapsed = time_elapsed+(t2-t1)
        
        # increase rho1 and rho2
        if np.mod(iter+2, rho_change_step)==0:
            rho1 = learning_fact*rho1
            rho2 = learning_fact*rho2
            gamma_val = max(gamma_val*gamma_factor,1)

        
        # evaluate this iteration
        temp1= (np.linalg.norm(x_sol-y1)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        temp2= (np.linalg.norm(x_sol-y2)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        if max(temp1,temp2)<=stop_threshold:
            print('iter: %d, stop_threshold: %.6f' %(iter, max(temp1,temp2)))
            break

        obj_list.append(compute_cost(x_sol,A,b))
        if len(obj_list)>=history_size:
            std_obj = compute_std_obj(obj_list, history_size)
            if std_obj<=std_threshold:
                print('iter: %d, std_threshold: %.6f' %(iter, std_obj))
                break
            
    
        cur_idx = x_sol>=0.5
        prev_idx = cur_idx
        cur_obj = compute_cost(prev_idx,A,b)
        
        # maintain best binary solution so far; in case the cost function oscillates
        if best_bin_obj >= cur_obj:
            best_bin_obj = cur_obj
            best_sol = x_sol
        
    return best_sol,x_sol,y1,y2,time_elapsed
    
    

"""
min_x x^T A x+b^Tx such that x is {0,1}^n; Cx=d
"""
def ADMM_bqp_linear_eq(A,b,C,d, all_params=None):

    initial_params = {'stop_threshold':1e-4,'std_threshold':1e-6,'gamma_val':1.6,'gamma_factor':0.95, 'rho_change_step':5, \
    'max_iters':1e3,'initial_rho':25,'history_size':3,'learning_fact':1+1/100,'x0':None,'pcg_tol':1e-4, 'pcg_maxiters':1e3,'rel_tol':5*1e-5, 'projection_lp':2}
    
    if all_params==None:
        all_params = initial_params
    else:
        for k in initial_params.keys():
            if k not in all_params.keys():
                all_params[k] = initial_params[k]
                            
    n = b.size
    stop_threshold = all_params['stop_threshold']
    std_threshold = all_params['std_threshold']
    max_iters = all_params['max_iters']
    initial_rho = all_params['initial_rho']
    rho_change_step = all_params['rho_change_step']
    gamma_val = all_params['gamma_val']
    learning_fact = all_params['learning_fact']
    history_size = all_params['history_size']
    projection_lp = all_params['projection_lp']
    gamma_factor = all_params['gamma_factor']
    pcg_tol = all_params['pcg_tol']
    pcg_maxiters = all_params['pcg_maxiters']


    # initialization
    x_sol = all_params['x0']
    y1 = x_sol    
    y2 = x_sol    
    z1 = np.zeros_like(y1)
    z2 = np.zeros_like(y2)
    z3 = np.zeros_like(d)
    rho1 = initial_rho
    rho2 = rho1
    rho3 = rho1
    obj_list = []
    Csq = C.transpose()@C  


    # initiate the binary solution
    prev_idx = x_sol
    best_sol = prev_idx
    best_bin_obj = compute_cost(best_sol,A,b)


    time_elapsed = 0
    for iter in range(int(max_iters)):
        t1 = time.time()
        
        # update y1: project onto box
        y1 = project_box(x_sol+z1/rho1)
    
    
        # update y2: project onto lp sphere
        y2 = project_shifted_Lp_ball(x_sol+z2/rho2, projection_lp)
        
     
        # update x: this is an exact solution to the subproblem
        # + solve a PD linear system, using pre-conditioned conjugate gradient algorithm  
        row = np.array(range(n))
        colum = np.array(range(n))
        data = (rho1+rho2)*np.ones(n)
        sparse_matrix = csc_matrix((data, (row, colum)), shape=(n, n))
        x_sol,cg_flag = linalg.cg(2*A+rho3*Csq+sparse_matrix, -(b+z1+z2+C.transpose()@z3)+rho1*y1+rho2*y2+rho3*C.transpose()@d, y1,  pcg_tol, pcg_maxiters)
        x_sol = x_sol.reshape(-1,1)
        
     
        # update z1 and z2 and z3
        z1 = z1+gamma_val*rho1*(x_sol-y1)
        z2 = z2+gamma_val*rho2*(x_sol-y2)
        z3 = z3+gamma_val*rho3*(C@x_sol-d)

        t2 = time.time()
        time_elapsed = time_elapsed+ (t2-t1)
        
        
        # increase rhos and update gamma is needed
        if np.mod(iter+1, rho_change_step)==0:
            rho1 = learning_fact*rho1
            rho2 = learning_fact*rho2
            rho3 = learning_fact*rho3
            gamma_val = max(gamma_val*gamma_factor,1)      
        
        
        # evaluate this iteration
        temp1= (np.linalg.norm(x_sol-y1)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        temp2= (np.linalg.norm(x_sol-y2)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        if max(temp1,temp2)<=stop_threshold:
            print('iter: %d, stop_threshold: %.6f' %(iter, max(temp1,temp2)))
            break

        obj_list.append(compute_cost(x_sol,A,b))
        if len(obj_list)>=history_size:
            std_obj = compute_std_obj(obj_list, history_size)
            if std_obj <= std_threshold:
                print('iter: %d, std_threshold: %.6f' %(iter, std_obj))
                break
        
        # maintain best binary solution so far; in case the cost function oscillates
        cur_idx = x_sol>=0.5
        prev_idx = cur_idx
        cur_obj = compute_cost(prev_idx,A,b)
    
        if best_bin_obj > cur_obj and max(temp1,temp2)>=1e-3:
            best_bin_obj = cur_obj
            best_sol = x_sol
            

    return best_sol,x_sol,y1,y2,time_elapsed



"""
This function implement lpbox ADMM to solve the following problem
min_x x'*A*x+b'*x such that x is {0,1}^n; Ex<=f
"""
def ADMM_bqp_linear_ineq(A,b,E,f,all_params):

    initial_params = {'stop_threshold':1e-4,'std_threshold':1e-6,'gamma_val':1.6,'gamma_factor':0.95, 'rho_change_step':5, \
        'max_iters':1e3, 'initial_rho':25, 'history_size':3, 'learning_fact':1+1/100, 'x0':[], 'pcg_tol':1e-4, 'pcg_maxiters':1e3, 'rel_tol':5e-5, 'projection_lp':2}

    if all_params==None:
        all_params = initial_params
    else:
        for k in initial_params.keys():
            if k not in all_params.keys():
                all_params[k] = initial_params[k]
                            
    n = b.size
    stop_threshold = all_params['stop_threshold']
    std_threshold = all_params['std_threshold']
    max_iters = all_params['max_iters']
    initial_rho = all_params['initial_rho']
    rho_change_step = all_params['rho_change_step']
    gamma_val = all_params['gamma_val']
    learning_fact = all_params['learning_fact']
    history_size = all_params['history_size']
    projection_lp = all_params['projection_lp']
    gamma_factor = all_params['gamma_factor']
    pcg_tol = all_params['pcg_tol']
    pcg_maxiters = all_params['pcg_maxiters']
    

    # initialization
    x_sol = all_params['x0']
    y1 = x_sol
    y2 = x_sol  
    y3 = f-E@x_sol
    z1 = np.zeros_like(y1)
    z2 = np.zeros_like(y2)
    z4 = np.zeros_like(y3)
    rho1 = initial_rho
    rho2 = rho1
    rho4 = rho1
    obj_list = []
    Esq = E.transpose()@E


    # initiate the binary solution
    prev_idx = x_sol
    best_sol = prev_idx
    best_bin_obj = compute_cost(best_sol,A,b)
    
    time_elapsed = 0;
    for iter in range(int(max_iters)):
        t1 = time.time()
        
        # update y1: project on box
        y1 = project_box(x_sol+z1/rho1)
        
        # update y2: project on circle
        y2 = project_shifted_Lp_ball(x_sol+z2/rho2, projection_lp)


        # update y3: project on non-negative quadrant
        y3 = f-E*x_sol-z4/rho4
        y3[y3<0]=0
        
     
        # update x: this is an exact solution to the subproblem
        # + solve a convex QP with linear constraints  
        
        row = np.array(range(n))
        colum = np.array(range(n))
        data = (rho1+rho2)*np.ones(n)
        sparse_matrix = csc_matrix((data, (row, colum)), shape=(n, n))
        x_sol,cg_flag = linalg.cg(2*A + rho4*Esq + sparse_matrix, \
                                  -(b+z1+z2+E.transpose()@z4)+rho1*y1+rho2*y2+rho4*E.transpose()@(f-y3), \
                                  x_sol, pcg_tol,pcg_maxiters)    
        x_sol = x_sol.reshape(-1,1)
        
    
        # update z1 and z2 and z4
        z1 = z1+gamma_val*rho1*(x_sol-y1)
        z2 = z2+gamma_val*rho2*(x_sol-y2)
        z4 = z4+gamma_val*rho4*(E*x_sol+y3-f)

        
        t2 = time.time()
        time_elapsed = time_elapsed+ (t2-t1)
        
        
        # increase rhos and update gamma is needed
        if np.mod(iter+2, rho_change_step)==0:
            rho1 = learning_fact*rho1
            rho2 = learning_fact*rho2
            rho4 = learning_fact*rho4
            gamma_val = max(gamma_val*gamma_factor,1)
       
        
        # evaluate this iteration
        temp1= (np.linalg.norm(x_sol-y1)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        temp2= (np.linalg.norm(x_sol-y2)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        if max(temp1,temp2)<=stop_threshold:
            print('iter: %d, stop_threshold: %.6f' %(iter, max(temp1,temp2)))
            break

        obj_list.append(compute_cost(x_sol,A,b))
        if len(obj_list)>=history_size:
            std_obj = compute_std_obj(obj_list, history_size)
            if std_obj <= std_threshold:
                print('iter: %d, std_threshold: %.6f' %(iter, std_obj))
                break
        
        # maintain best binary solution so far; in case the cost function oscillates
        cur_idx = x_sol>=0.5
        prev_idx = cur_idx
        cur_obj = compute_cost(prev_idx,A,b)
    
        if best_bin_obj > cur_obj and np.linalg.norm(E@x_sol-f)<=((1e-3)*np.sqrt(n)):
            best_bin_obj = cur_obj
            best_sol = x_sol
        
        
    return best_sol,x_sol,y1,y2,time_elapsed


"""
This function implement lpbox ADMM to solve the following problem
min_x x'*A*x+b'*x such that x is {0,1}^n; Cx=d, Ex<=f
"""
def ADMM_bqp_linear_eq_and_uneq(A,b,C,d,E,f, all_params=None):

    initial_params = {'stop_threshold':1e-4,'std_threshold':1e-6,'gamma_val':1.6,'gamma_factor':0.95, 'rho_change_step':5, \
    'max_iters':1e3,'initial_rho':25,'history_size':3,'learning_fact':1+1/100,'x0':None,'pcg_tol':1e-4, 'pcg_maxiters':1e3,'rel_tol':5*1e-5, 'projection_lp':2}
    
    if all_params==None:
        all_params = initial_params
    else:
        for k in initial_params.keys():
            if k not in all_params.keys():
                all_params[k] = initial_params[k]
                            
    n = b.size
    stop_threshold = all_params['stop_threshold']
    std_threshold = all_params['std_threshold']
    max_iters = all_params['max_iters']
    initial_rho = all_params['initial_rho']
    rho_change_step = all_params['rho_change_step']
    gamma_val = all_params['gamma_val']
    learning_fact = all_params['learning_fact']
    history_size = all_params['history_size']
    projection_lp = all_params['projection_lp']
    gamma_factor = all_params['gamma_factor']
    pcg_tol = all_params['pcg_tol']
    pcg_maxiters = all_params['pcg_maxiters']


    # initialization
    x_sol = all_params['x0']
    y1 = x_sol    
    y2 = x_sol
    y3 = f-E@x_sol   
    
    z1 = np.zeros_like(y1)
    z2 = np.zeros_like(y2)
    z3 = np.zeros_like(d)
    z4 = np.zeros_like(y3)
    
    rho1 = initial_rho
    rho2 = rho1
    rho3 = rho1
    rho4 = rho1
    
    obj_list = []
    Csq = C.transpose()@C  
    Esq = E.transpose()@E


    # initiate the binary solution
    prev_idx = x_sol
    best_sol = prev_idx
    best_bin_obj = compute_cost(best_sol,A,b)


    time_elapsed = 0
    for iter in range(int(max_iters)):
        t1 = time.time()
        
        # update y1: project onto box
        y1 = project_box(x_sol+z1/rho1)
    
    
        # update y2: project onto lp sphere
        y2 = project_shifted_Lp_ball(x_sol+z2/rho2, projection_lp)
        
        
        # update y3: project on non-negative quadrant
        y3 = f-E@x_sol-z4/rho4
        y3[y3<0]=0
        
        # update x: this is an exact solution to the subproblem
        # + solve a PD linear system, using pre-conditioned conjugate gradient algorithm  
        row = np.array(range(n))
        colum = np.array(range(n))
        data = (rho1+rho2)*np.ones(n)
        sparse_matrix = csc_matrix((data, (row, colum)), shape=(n, n))
        

        x_sol,cg_flag = linalg.cg(2*A+rho3*Csq+rho4*Esq+sparse_matrix, -(b+z1+z2+C.transpose()@z3+E.transpose()@z4) \
                                                                       +rho1*y1+rho2*y2+rho3*C.transpose()@d+rho4*E.transpose()@(f-y3), \
                                                                         y1,  pcg_tol, pcg_maxiters)
        x_sol = x_sol.reshape(-1,1)
        
     
        # update z1 and z2 and z3 and z4
        z1 = z1+gamma_val*rho1*(x_sol-y1)
        z2 = z2+gamma_val*rho2*(x_sol-y2)
        z3 = z3+gamma_val*rho3*(C@x_sol-d)
        z4 = z4+gamma_val*rho4*(E@x_sol+y3-f)

        t2 = time.time()
        time_elapsed = time_elapsed+ (t2-t1)
        
        
        # increase rhos and update gamma is needed
        if np.mod(iter+1, rho_change_step)==0:
            rho1 = learning_fact*rho1
            rho2 = learning_fact*rho2
            rho3 = learning_fact*rho3
            rho4 = learning_fact*rho4
            gamma_val = max(gamma_val*gamma_factor,1)      
        
        
        # evaluate this iteration
        temp1= (np.linalg.norm(x_sol-y1)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        temp2= (np.linalg.norm(x_sol-y2)) / max(np.linalg.norm(x_sol), 2.2204e-16)
        if max(temp1,temp2)<=stop_threshold:
            print('iter: %d, stop_threshold: %.6f' %(iter, max(temp1,temp2)))
            break

        obj_list.append(compute_cost(x_sol,A,b))
        if len(obj_list)>=history_size:
            std_obj = compute_std_obj(obj_list, history_size)
            if std_obj <= std_threshold:
                print('iter: %d, std_threshold: %.6f' %(iter, std_obj))
                break
        
        # maintain best binary solution so far; in case the cost function oscillates
        cur_idx = x_sol>=0.5
        prev_idx = cur_idx
        cur_obj = compute_cost(prev_idx,A,b)
    
        if best_bin_obj > cur_obj and max(temp1,temp2)>=1e-3:
            best_bin_obj = cur_obj
            best_sol = x_sol
            

    return best_sol,x_sol,y1,y2,time_elapsed



def project_box(x):
    xp = x
    xp[x>1]=1
    xp[x<0]=0

    return xp


def project_shifted_Lp_ball(x, p):
    shift_vec = 1/2*np.ones((x.size, 1))
    shift_x = x-shift_vec
    normp_shift = np.linalg.norm(shift_x, p)
    n = x.size
    xp = (n**(1/p)) * shift_x / (2*normp_shift) + shift_vec

    return xp

    
def compute_cost(x,A,b):
    c = x.transpose()@A@x + b.transpose()@x
    
    return c


def compute_std_obj(obj_list, history_size):
    
    std_obj = np.std(obj_list[-1-history_size:])
    
    # normalize 
    std_obj = std_obj/abs(obj_list[-1])

        
    return std_obj[0][0]
