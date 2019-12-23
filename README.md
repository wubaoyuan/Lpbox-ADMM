# Lp-Box ADMM for Integer Programming

This project provides two implementations of [Lp-Box ADMM](https://ieeexplore.ieee.org/document/8378001/) using
* [Matlab](matlab): full codes and full demos to reproduce all reported results in the manuscript.
* [Python](python): full codes and one simple demo to demonstrate the usage. 

## Basic idea

[Lpbox-ADMM](https://ieeexplore.ieee.org/document/8378001/) is a general optimization method used for interger programming. Any problems with discrete constraints can be solved using Lp-Box ADMM.  Since any discrete constraint can be easily transformed to binary constraint with an additional simplex constraint, we focus on the following problem with binary constraints:
$$
  \mathop{\min}_x \ f(x) \quad \text{s.t.} \quad x \in \{0,1\}^n, x \in \mathcal{C}
$$

To solve the binary discrete constraint, we propose to replace it with an equivalent set of continuous constraints.  
$$
  x \in \{0,1\}^n \leftrightarrow x\in\[0,1\]^n \cap \{||x-\frac{1}{2}||_p^p=\frac{n}{2^p}\}
$$  

where $\[0,1\]^n$ is called box-constraint, and $\{||x-\frac{1}{2}||_p^p=\frac{n}{2^p}\}$ is denoted as $\ell_p$-norm constraint.
The geometric illustration of the equivalence between lpbox intersection and the set of binary points is presented as follows. For clarity, we just show the cases when p={1,2,5}

<div align="center">
<img src="/lpbox.png">
</div>


## Binary Quadratic Programming (BQP) 
In this project, we specify the optimization object as the QP problem with binary discrete constraint. The BQP problem can be denoted as follows:
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1, C_2 x \leq d_2
$$

which includes binary discrete constraint, linear equation constraint and linear inequation constraint. 

## Python usages
To briefly show how to use this algorithm, we also present three python interfaces for the BQP problems.  
1. unconstrained BQP  
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n
$$  

2. BQP with linear equality constraints 
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1
$$

3. BQP with linear inequality constraints 
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_2 x \leq d_2
$$  

4. BQP with linear equality and inequality constraints 
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1, C_2 x \leq d_2
$$  

To use these three interfaces, you need to firstly add these codes to you project.
```
from functions.lpbox_admm import ADMM_bqp_unconstrained
from functions.lpbox_admm import ADMM_bqp_linear_eq
from functions.lpbox_admm import ADMM_bqp_linear_ineq
from functions.lpbox_admm import ADMM_bqp_linear_eq_and_ineq
```

## Demo
To help you quickly understand the how to use this algorithm, we present a simple demo of image segmentation by solving unconstrained BQP problem. The demo calls the ```ADMM_bqp_unconstrained``` function to solve the unconstrained BQP problem, you can run the following code to have a try.
```
python demo_image_segmentation.py
```
We randomly initalize the image with 0,1 and optimize the image to the segmentation result.

<div align="center">
<img src="/demo/show_image.png">
</div>

## Other applications
1. #### [Model compression](https://github.com/wubaoyuan/CNN-FCF-CVPR-2019)

This work proposes to conduct filter selection and filter learning in a unified model. They define a factorized convolutional filter (FCF), consisting of a standard real-valued convolutional filter and a binary scalar, as well as a dot-product operator between them. They train a CNN model with factorized convolutional filters (CNN-FCF) by updating the standard filter using back-propagation, while updating the binary scalar using the alternating direction method of multipliers (ADMM) based optimization method.

## Citation
If you adopt the Lp-Box ADMM algorithm in your project, please cite as follows.
```
@article{wu2018lp,
  title={lp-box ADMM: A versatile framework for integer programming},
  author={Wu, Baoyuan and Ghanem, Bernard},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2018},
  publisher={IEEE}
}
```
If any problem, feel free to contact with [Baoyuan Wu](https://sites.google.com/site/baoyuanwu2015/home) (wubaoyuan1987@gmail.com).
