# Lp-Box ADMM for Integer Programming

[Lpbox-ADMM](https://ieeexplore.ieee.org/document/8378001/) is a general optimization method for ANY interger programming (IP), and it has been published in TPAMI 2018.  
This project provides two implementations of Lp-Box ADMM:
* [Matlab](matlab): full codes and full demos to reproduce all reported results in the manuscript.
* [Python](python): full codes and one simple demo to demonstrate the usage. 

**IMPORTANT NOTE:** We have just provided the implementation of Lp-Box ADMM for the BQP problem. However, Lp-Box ADMM can be naturally applied to ANY IP tasks. Welcome to contribute your IP tasks and implementations to this project. 

## Contents
  
* [Basic idea](#basic-idea)
  * [Equivalent replacement of binary constraint](#equivalent-replacement-of-binary-constraint)
  * [Constraint splitting via extra variables](#constraint-splitting-via-extra-variables)
  
* [Python usage](#python-usage)
  * [Binary quadratic programming](#binary-quadratic-programming)
  * [Demo of image segmentation](#demo-of-image-segmentation)
  
* [Applications and extensions](#applications-and-extensions)
  * [1 Deep model compression](#1-deep-model-compression)
  * [2 MAP inference for probabilistic graphical models](#2-map-inference-for-probabilistic-graphical-models)
  * [3 Kmeans clustering](#3-kmeans-clustering)
  * [4 Others](#4-others)
  
* [Citation](#citation)


## [Basic idea](#basic-idea)
[[back to top](#)]

Since any discrete constraint can be easily transformed to the binary constraint with an additional simplex constraint, in the following we focus on the following problem with binary constraints:
$$
  \mathop{\min}_x \ f(x) \quad \text{s.t.} \quad x \in \{0,1\}^n, x \in \mathcal{C}
$$

#### [Equivalent replacement of binary constraint](#equivalent-replacement-of-binary-constraint)
We propose to replace the binary constraint with an equivalent set of continuous constraints.  
$$
  x \in \{0,1\}^n \leftrightarrow x\in\[0,1\]^n \cap \{||x-\frac{1}{2}||_p^p=\frac{n}{2^p}\}
$$  

where $\mathcal{S}_b = \[0,1\]^n$ is called box-constraint, and $\mathcal{S}_p = \{||x-\frac{1}{2}||_p^p=\frac{n}{2^p}\}$ is denoted as $\ell_p$-norm constraint.

The geometric illustration of the constraint equivalence is presented as follows. For clarity, we just show the cases when p={1,2,5}.
<div align="center">
<img src="/lpbox.png">
</div>

#### [Constraint splitting via extra variables](#constraint-splitting-via-extra-variables)
We further introduce two extra variables $y_1$ and $y_2$ to split the constraints onto different variables, such that the constraints can be iteratively/gradually satisfied, as follows
$$
  \mathop{\min}_x \ f(x) \quad \text{s.t.} \quad x \in \{0,1\}^n, x \in \mathcal{C}, x=y_1, x=y_2, y_1 \in \mathcal{S}_b, y_2 \in \mathcal{S}_p.
$$

The above problem can be easily solved by the alternating direction method of multipliers (ADMM). 


## [Python usage](#python-usage)
[[back to top](#)]

### [Binary quadratic programming](#binary-quadratic-programming)
Since many important applications can be formulated as BQP, in this project we present the demo of using Lp-Box ADMM to solve the BQP problem, which is formulated as follows
$$
  \mathop{\min}_x \ x^\top Ax+b^\top \quad \text{s.t.} \ x \in \{0,1\}^n, C_1 x=d_1, C_2 x \leq d_2
$$

which includes the binary constraint, linear equality and inequality constraint. 
To facilitate the usages, we provide four python functions for the BQP problems with different constraints.  
To use these functions, you just need to import them to your demo, as follows
```
from functions.lpbox_admm import ADMM_bqp_unconstrained
from functions.lpbox_admm import ADMM_bqp_linear_eq
from functions.lpbox_admm import ADMM_bqp_linear_ineq
from functions.lpbox_admm import ADMM_bqp_linear_eq_and_ineq
```

<!---
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
-->


### [Demo of image segmentation](#demo-of-image-segmentation)
We present a simple demo of image segmentation by solving unconstrained BQP problem, which calls the ```ADMM_bqp_unconstrained``` function, as follows
```
python demo_image_segmentation.py
```
The randomly initialized image and the segmentation result are shown as follows
<div align="center">
<img src="/python/demo/show_image.png">
</div>

## [Applications and extensions](#applications-and-extensions)
[[back to top](#)]

#### [1 Deep model compression](#1-deep-model-compression)

This work has been published in CVPR 2019, "Compressing Convolutional Neural Networks via Factorized Convolutional Filters" ([pdf](http://openaccess.thecvf.com/content_CVPR_2019/papers/Li_Compressing_Convolutional_Neural_Networks_via_Factorized_Convolutional_Filters_CVPR_2019_paper.pdf), [github](https://github.com/wubaoyuan/CNN-FCF-CVPR-2019))

We applied the idea of Lp-Box ADMM to deep model compression, which learns and selects the convolutional filters in a unified model. Specifically, we fitstly define a factorized convolutional filter (FCF), consisting of a standard real-valued convolutional filter and a binary selection scalar, as well as a dot-product operator between them. Then, we train CNN model with factorized convolutional filters (CNN-FCF), by updating the standard filter using back-propagation, while updating the binary scalar using the alternating direction method of multipliers (ADMM) based optimization method. The framework of the standard filter pruning (top) and the proposed CNN-FCF based pruning (bottom) are shown in the following figure.

<div align="center">
<img src="/figures/CNN-FCF.png">
</div>


#### [2 MAP inference for probabilistic graphical models](#2-map-inference-for-probabilistic-graphical-models)

This work has been accepted to IJCV, "MAP Inference via L2-Sphere Linear Program Reformulation" ([Arxiv](https://arxiv.org/pdf/1905.03433.pdf), the github project will be released soon).

MAP inference is a fundamental task in probabilistic graphical models, which aims to infer the most probable label configuration of a probabilistic graphical model (e.g., MRF, CRF, HMM). MAP inference can be formulated as an integer programming, based on the factor graph (any graphical model can be transformed to a corresponding factor graph), as follows
$$
  \text{MAP}(\boldsymbol{\theta}) = \text{ILP}(\boldsymbol{\theta}) = \mathop{\max}_{\boldsymbol{\mu}, \boldsymbol{v}} < \boldsymbol{\theta}, \boldsymbol{\mu} > 
  ~  \quad \text{s.t.} \quad \boldsymbol{\mu} \in \mathcal{L}_G \cap \{0, 1\}^{|\boldsymbol{\mu}|}.
$$

<!---
* MAP inference of a MRF model
\begin{equation}
  \mathop{\min}_x  \text{MAP}(\boldsymbol{\theta}) = \max_{\mathbf{G} \in \mathcal{X}} \sum_{i\in V} \boldsymbol{\theta}_i(\mathbf{g}_i) + \sum_{\alpha \in F} \boldsymbol{\theta}_{\alpha}(\mathbf{g}_{\alpha})
\end{equation}
* MAP inference as integer programming 
$$
  \mathop{\min}_x \text{ILP}(\boldsymbol{\theta}) = \max_{\boldsymbol{\mu}} \sum_{i\in V} \boldsymbol{\theta}_i^\top \boldsymbol{\mu}_i + \sum_{\alpha \in F} \boldsymbol{\theta}_{\alpha}^\top \boldsymbol{\mu}_{\alpha}
 = \max_{\boldsymbol{\mu}}  \langle \boldsymbol{\theta}, \boldsymbol{\mu} \rangle, ~
 \text{s.t.} ~  \boldsymbol{\mu} \in \mathcal{L}_G \cap \{0,1\}^{|\boldsymbol{\mu}|},
$$ 
$$
  \mathop{\min}_x \ \text{ILP}(\boldsymbol{\theta}) = \max_{\boldsymbol{\mu}} 
$$  
--->

Inspired by the idea of Lp-Box ADMM, we firstly remove the binary constraint, while adding the $\ell_2$-sphere constraint onto the variable or the factor nodes. Then, we introduce an extra variable/factor nodes to split the $\ell_2$-sphere constraint. 
$$
  \text{LS-LP}(\boldsymbol{\theta}) = \mathop{\max}_{\boldsymbol{\mu}, \boldsymbol{v}} < \boldsymbol{\theta}, \boldsymbol{\mu} > 
  ~  \quad \text{s.t.} \quad \boldsymbol{\mu} \in \mathcal{L}_G, \boldsymbol{v} \in \mathcal{S}, \boldsymbol{\mu}_i = \boldsymbol{v}_i, 
  i \in V.
$$

It is easily proved that it is an equivalent continuous reformulation to the original integer programming of MAP inference, which was then efficiently solved by ADMM. And, we prove that it is globally convergent to epsilon-KKT solution.

#### [3 Kmeans clustering](#3-kmeans-clustering)

#### [4 Others](#4-others)
      

## [Citation](#citation)
[[back to top](#)]

If any problem, feel free to contact with [Baoyuan Wu](https://sites.google.com/site/baoyuanwu2015/home) (wubaoyuan1987@gmail.com).

```
@article{wu2018lp,
  title={lp-box ADMM: A versatile framework for integer programming},
  author={Wu, Baoyuan and Ghanem, Bernard},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2018},
  publisher={IEEE}
}
```

