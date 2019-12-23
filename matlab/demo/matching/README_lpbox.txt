
============
To evaluate our Lp-box ADMM algorithm on graph matching task, we utilize the toolbox developed by Feng Zhou. 
Please check the demo file 'demo_graph_matching_LpADMM_and_other_methods.m'
The original readme file written by Feng Zhou is as follows. 

--- Baoyuan Wu (wubaoyuan1987@gmail.com), 2018/06/08
============



---------------- The following is the original README written by Feng Zhou -----------------

============
Introduction
============
    This page contains software and instructions for factorized graph matching (FGM) [1] [2].
    In addition, we include the following well-known methods as baselines:

	  spectral matching (SM) [3],
      spectral matching with affine constraints (SMAC) [4],
      graduated assignment (GA) [5],
	  probabilistic matching (PM) [6],
	  integer projected fixed point method (IPFP) [7],
      re-weighted random walk matching (RRWM) [8].

    The implementations of the above methods are taken from the authors' websites.


============
Installation
============
    1. unzip "fgm.zip" to your folder;
    2. Run "make" to compile all C++ files;
    3. Run "addPath" to add sub-directories into the path of Matlab.
    4. Run "demoXXX" or "testXXX".


============
Instructions
============
    The package of "fgm.zip" contains the following files and folders:

    ./data: This folder contains the CMU House Image dataset.

    ./save: This folder contains the experimental results reported in the paper.

    ./src: This folder contains the main implmentation of FGM as well as other baselines.

    ./lib: This folder contains some necessary library functions.

    ./make.m: Matlab makefile for C++ code.

    ./addPath.m: Adds the sub-directories into the path of Matlab.

    ./demoToy.m: A demo comparison of different graph matching methods on the synthetic dataset.

    ./demoHouse.m: A demo comparison of different graph matching methods on the on CMU House dataset.
  
    ./testToy.m: Testing the performance of different graph matching methods on the synthetic dataset.
                 This is a similar function used for reporting (Fig. 4) the first experiment (Sec 5.1) in the CVPR 2012 paper

    ./testHouse.m Testing the performance of different graph matching methods on CMU House Image dataset.
			      This is the same function used for reporting (Fig. 4) the first experiment (Sec 5.1) in the CVPR 2013 paper.


========
C++ Code
========
    In order to achieve efficiency, we provide several C++ codes under
    "src/asg/fgm/matrix" for performing matrix products between
    binary matrices. For instance, the function "multiGXH.cpp" is used
    to more efficiently compute the matrix product, "G^T * X * H", where G and H
    are two binary matrices.


==========
References
==========
    [1] F. Zhou and F. De la Torre, "Deformable Graph Matching," in IEEE
    Conference on Computer Vision and Pattern Recognition (CVPR), 2013.
    
    [2] F. Zhou and F. De la Torre, "Factorized Graph Matching," in IEEE
    Conference on Computer Vision and Pattern Recognition (CVPR), 2012.
    
    [3] M. Leordeanu and M. Hebert, "A spectral technique for
    correspondence problems using pairwise constraints," in International
    Conference on Computer Vision (ICCV), 2005.
    
    [4] T. Cour, P. Srinivasan and J. Shi, "Balanced Graph Matching", in
    Advances in Neural Information Processing Systems (NIPS), 2006
    
    [5] S. Gold and A. Rangarajan, "A Graduated Assignment Algorithm for
    Graph Matching", IEEE Transactions on Pattern Analysis and Machine
    Intelligence (PAMI), 1996
    
    [6] R. Zass and A. Shashua, "Probabilistic Graph and Hypergraph
    Matching", in IEEE Conference on Computer Vision and Pattern
    Recognition (CVPR), 2008
    
    [7] M. Leordeanu, M. Hebert and R. Sukthankar, "An Integer Projected
    Fixed Point Method for Graph Matching and MAP Inference", in Advances
    in Neural Information Processing Systems (NIPS), 2009
    
    [8] M. Cho, J. Lee and K. Lee, "Reweighted Random Walks for Graph
	Matching", in European Conference on Computer Vision (ECCV), 2010









Copyright
    This software is free for use in research projects. If you
    publish results obtained using this software, please use this
    citation.
@inproceedings{ZhouD12b,
   author   = {Feng Zhou and Fernando {De la Torre}},
   title    = {Factorized Graph Matching},
   booktitle    = {IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
   year     = {2012},
}
If you have any question, please feel free to contact Feng Zhou (zhfe99@gmail.com).
