# Lp-Box ADMM for Integer Programming

This project provides two implementations of [Lp-Box ADMM](https://ieeexplore.ieee.org/document/8378001/) using
* [Matlab](matlab): full codes and full demos to reproduce all reported results in the manuscript.
* [Python](python): full codes and one simple demo to demonstrate the usage. 

## Basic idea

[Lpbox-ADMM](https://ieeexplore.ieee.org/document/8378001/) is a general optimization method for ANY interger programming. Since any discrete constraint can be easily transformed to the binary constraint with an additional simplex constraint, we focus on the following problem with binary constraints:
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


## Binary quadratic programming (BQP) 
In this project, we specify the optimization object as the QP problem with binary constraints. The BQP problem can be denoted as follows:
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1, C_2 x \leq d_2
$$

which includes the binary constraint, linear equality and inequality constraint. 

## Python usages

### Functions
To facilitate the usages of the Lp-Box ADMM method, we provide four python functions for the BQP problems with different constraints.  
* unconstrained BQP  
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n
$$  


* BQP with linear equality constraints 
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1
$$

* BQP with linear inequality constraints 
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_2 x \leq d_2
$$  

* BQP with linear equality and inequality constraints 
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1, C_2 x \leq d_2
$$  

To use these functions, you need to import them to your demo, as follows
```
from functions.lpbox_admm import ADMM_bqp_unconstrained
from functions.lpbox_admm import ADMM_bqp_linear_eq
from functions.lpbox_admm import ADMM_bqp_linear_ineq
from functions.lpbox_admm import ADMM_bqp_linear_eq_and_ineq
```

### Demo
We present a simple demo of image segmentation by solving unconstrained BQP problem, which calls the ```ADMM_bqp_unconstrained``` function, as follows
```
python demo_image_segmentation.py
```
The randomly initialized image and the segmentation result are shown as follows
<div align="center">
<img src="/python/demo/show_image.png">
</div>

## Aplications and extensions
1. #### [Model compression for deep neural networks](https://github.com/wubaoyuan/CNN-FCF-CVPR-2019)

We applied the idea of Lp-Box ADMM to deep model compression, which learns and selects the convolutional filters in a unified model. Specifically, we fitstly define a factorized convolutional filter (FCF), consisting of a standard real-valued convolutional filter and a binary selection scalar, as well as a dot-product operator between them. Then, we train CNN model with factorized convolutional filters (CNN-FCF), by updating the standard filter using back-propagation, while updating the binary scalar using the alternating direction method of multipliers (ADMM) based optimization method. The framework of the standard filter pruning (top) and the proposed CNN-FCF based pruning (bottom) are shown in the following figure.

<div align="center">
<img src="/figures/CNN-FCF.png">
</div>

This work has been published in CVPR 2019, "Compressing Convolutional Neural Networks via Factorized Convolutional Filters" ([pdf](http://openaccess.thecvf.com/content_CVPR_2019/papers/Li_Compressing_Convolutional_Neural_Networks_via_Factorized_Convolutional_Filters_CVPR_2019_paper.pdf), [github](https://github.com/wubaoyuan/CNN-FCF-CVPR-2019))


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
