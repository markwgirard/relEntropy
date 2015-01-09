relEntropy
==========
This is a MATLAB routine for computing the relative entropy of entanglement of input bipartite density matrices.

Standard usage:

    [Xopt,relEnt] = relEntropy(m,n,rho);
where 
* `m`,`n` are dimensions of subspaces
* `rho` is input `mn x mn` density matrix
* `Xopt` is closest PPT state to `rho`
* `relEnt` is numerical approximation to the relative entropy of entanglement of `rho`

## Requirements

* MATLAB
* CVX by CVX Research, Inc. (http://cvxr.com/cvx/)
* PartialTranspose routine from QETLAB (http://www.qetlab.com) (see: (https://github.com/nathanieljohnston/QETLAB/blob/master/PartialTranspose.m))

## References
Y. Zinchenko, S. Friedland, and G. Gour. *Numerical estimation of the relative entropy of entanglement*.
[Phys. Rev. A 82, 052336](http://dx.doi.org/10.1103/PhysRevA.82.052336) (2010). 


