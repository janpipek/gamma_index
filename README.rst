gamma_index
===========

Small Python library for calculation of gamma index.

The module offers two functions:

`gamma_matrix` 
    This function calculates gamma index in each point for two
    n-dimensional distributions (with a same shape). It is not
    very optimized, it's O(*n* x *n*), where *n* is the total
    number of elements of the matrix. It takes minutes to calculate
    the 3-D matrix for 50x50x50 distribution. Time scales
    roughly as n^6 (n = size of square matrix in one dimension)

`gamma_matrix_c`
	This functiion calculates gamma index using a slightly optimized
	C algorithm. It becomes more efficient than gamma_matrix
	with a matrix size of roughly 1000 elements for default parameters.
	However, the relative efficiency of this function highly depends
	on the parameters and values (the better agreement, the faster the algorithm is).
	Time scales roughly as n^3.2 (n = size of square matrix in one dimension)

`gamma_matrix_pass`
    This function returns a boolean matrix (passed / not passed) 
    with optimized search only on a narrow neighbourhood of
    each point. 


TODO
----

* Prepare an optimized version of `gamma_matrix` based on 
  layers
* Add tests