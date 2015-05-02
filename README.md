# gamma_index

Small Python library for calculation of gamma index.

The module offers two functions:

`gamma_matrix` 
: This function calculates gamma index in each point for two
  n-dimensional distributions (with a same shape). It is not
  very optimized, it's O(*n* x *n*), where *n* is the total
  number of elements of the matrix. It takes minutes to calculate
  the 3-D matrix for 50x50x50 distribution.

`gamma_matrix_pass`
: This function returns a boolean matrix (passed / not passed) 
  with optimized search only on a narrow neighbourhood of
  each point. 

## TODO

* Prepare an optimized version of `gamma_matrix` based on 
  layers
* Add tests