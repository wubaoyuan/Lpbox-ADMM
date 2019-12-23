# Lp-Box ADMM for Integer Programming

Lpbox-ADMM ([main manuscript](https://ieeexplore.ieee.org/document/8378001/), [supplementary material](https://www.dropbox.com/s/er6nwiia1c9734t/Lpbox_ADMM_PAMI_revision_supplementary.pdf?dl=0))  is a generic optimization method for interger programming (IP), and it has been published in TPAMI 2018. This project provides two implementations of Lp-Box ADMM:
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
<img src="/lpbox.png"  width="500"/>
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
<img src="/python/demo/show_image.png"  width="600"/>
</div>

## [Applications and extensions](#applications-and-extensions)
[[back to top](#)]

#### [1 Deep model compression](#1-deep-model-compression)

This work has been published in CVPR 2019, "Compressing Convolutional Neural Networks via Factorized Convolutional Filters" ([pdf](http://openaccess.thecvf.com/content_CVPR_2019/papers/Li_Compressing_Convolutional_Neural_Networks_via_Factorized_Convolutional_Filters_CVPR_2019_paper.pdf), [github](https://github.com/wubaoyuan/CNN-FCF-CVPR-2019))

We applied the idea of Lp-Box ADMM to deep model compression, which learns and selects the convolutional filters in a unified model. Specifically, we fitstly define a factorized convolutional filter (FCF), consisting of a standard real-valued convolutional filter and a binary selection scalar, as well as a dot-product operator between them. Then, we train CNN model with factorized convolutional filters (CNN-FCF), by updating the standard filter using back-propagation, while updating the binary scalar using the alternating direction method of multipliers (ADMM) based optimization method. The framework of the standard filter pruning (top) and the proposed CNN-FCF based pruning (bottom) are shown in the following figure.

<div align="center">
<img src="/figures/CNN-FCF.png"  width="600"/>
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

It is easily proved that $ \text{LS-LP}(\boldsymbol{\theta}) = \text{MAP}(\boldsymbol{\theta}) = \text{ILP}(\boldsymbol{\theta})$.  LS-LP can be efficiently solved by ADMM, which is proved to be globally convergent to epsilon-KKT solution of the original MAP inference.



<div align="center">
<img src="/figures/factor-graph.png"  width="500"/>
</div>


#### [3 Kmeans clustering](#3-kmeans-clustering)

This work is presented in Arxiv, "Constrained K-means with General Pairwise and Cardinality Constraints"([Arxiv](https://arxiv.org/pdf/1907.10410.pdf)). 

K-means is one of the most popular classic clustering algorithms. However, the orginal K-means is unstable. One enhanced approach is inserting some user preferences (e.g., pairwise constraints) into K-means, using some heuristic strategies. [One recent work](http://www.optimization-online.org/DB_FILE/2005/04/1114.pdf) formulates K-means as an integer programming. Based on this formulation, different types of user preferences can be naturally embedded as constraints, such as cardinality constraints, must/cannot-link constraints. We adopt the Lp-Box ADMM algorithm to optimize this IP problem. 


#### [4 Others](#4-others)

The Lp-Box ADMM has been applied to many types of applications, such as [hash code learning](http://cfm.uestc.edu.cn/~fshen/SADH.pdf), 
[low-density parity-check (LDPC)](https://arxiv.org/pdf/1711.10767.pdf), [feature selection](https://www.ijcai.org/proceedings/2017/0228.pdf), [data hiding](http://www.busim.ee.boun.edu.tr/~sankur/SankurFolder/Conf_EUSIPCO_2018_Robust%20Data%20Hiding.pdf), etc.

Looking forward to more applications using Lp-Box ADMM. 

      

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

