__all__ = [ "gamma_matrix", "gamma_matrix_pass"]


import numpy as np


def gamma_matrix(rm, tm, dta=1.0, dd=0.05):
    '''Compute matrix of gammma indices.

    :param rm: reference matrix (relative values assumed)
    :param tm: tested matrix (relative values assumed)
    :param dta: maximum distance-to-agreement (in voxels)
    :param dd: maximum dose difference (absolute, not percent!)

    :type rm: numpy.ndarray
    :type tm: numpy.ndarray
    :type dta: float
    :type dd: float
    :rtype: numpy.ndarray

    It can be evaluated on matrices of any dimensionality.
    '''
    # Check validity of input
    if rm.shape != tm.shape:
        raise Exception("Cannot compute for matrices of different sizes.")

    # Result matrix
    output = np.ndarray(rm.shape, dtype=np.float64)

    # Help scaling variables
    dta_scale = dta ** 2
    dd_scale = dd ** 2

    # Index matrices
    indices = np.indices(rm.shape)

    it = np.nditer(rm, ("multi_index",))
    while not it.finished:
        index = it.multi_index

        # Dose difference to every point (squared)
        dd2 = (tm - it.value) ** 2

        # Distance to every point (squared)
        dist2 = np.sum((indices - np.array(index).reshape(len(rm.shape),1,1,1)) ** 2, axis=0)

        # Minimum of the sum
        output[index] = np.sqrt(np.nanmin(dd2 / dd_scale + dist2 / dta_scale))
        it.iternext()
    return output


def gamma_matrix_pass(rm, tm, dta=1.0, dd=0.05, ignore=lambda value: False):
    '''Effectively check which points in a matrix pair pass gamma index test.

    :param rm: reference matrix (relative values assumed)
    :param tm: tested matrix (relative values assumed)
    :param dta: maximum distance-to-agreement (in voxels)
    :param dd: maximum dose difference
    :param ignore: function called on dose in rm values
        if it returns True, point is ignored (gammma <- np.nan)
    :rtype: numpy.ndarray

    It can be evaluated on matrices of any dimensionality.
    Optimized in that only surrounding region that can possibly
        pass dta criterion is checked for each point.
    '''
    # Check validity of input
    if rm.shape != tm.shape:
        raise Exception("Cannot compute for matrices of different sizes.")

    shape = rm.shape
    ndim = rm.ndim

    # Result matrix
    output = np.ndarray(rm.shape, dtype=bool)

    # Help scaling variables
    dta_scale = dta ** 2
    dd_scale = dd ** 2

    # Index matrices
    indices = np.indices(rm.shape)

    # How many points (*2 + 1)
    npoints = int(dta)

    it = np.nditer(rm, ("multi_index",))
    while not it.finished:
        index = tuple(it.multi_index)

        if ignore(it.value):
            output[index] = np.nan
            it.iternext()
            continue

        slices = [ slice(max(0, index[i] - npoints), min(shape[i], index[i] + npoints + 1)) for i in range(ndim) ]
        subtm = tm[slices]

        # Dose difference to every point (squared)
        dd2 = (subtm - it.value) ** 2

        # Distance to every point (squared)
        dist2 = np.sum((indices[[slice(None, None)] + slices] - np.array(index).reshape(ndim,1,1,1)) ** 2, axis=0)

        # Minimum of the sum
        output[index] = np.sqrt(np.nanmin(dd2 / dd_scale + dist2 / dta_scale)) < 1.0
        it.iternext()
    return output


def gamma_matrix_c(rm, tm, dta=1.0, dd=0.05):
    '''Compute matrix of gammma indices (faster method using C)

    :param rm: reference matrix (relative values assumed)
    :param tm: tested matrix (relative values assumed)
    :param dta: maximum distance-to-agreement (in voxels)
    :param dd: maximum dose difference (absolute, not percent!)

    :type rm: numpy.ndarray
    :type tm: numpy.ndarray
    :type dta: float
    :type dd: float
    :rtype: numpy.ndarray

    It can be evaluated on matrices of any dimensionality (but of a same shape).
    Requires previously built libgamma.so (see Makefile) and works only in Linux.
    '''
    import ctypes
    import os

    so_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "libgamma.so")

    dll = ctypes.CDLL(so_path)
    gamma_index_f = dll.gamma_index

    # Check validity of input
    if rm.shape != tm.shape:
        raise Exception("Cannot compute for matrices of different sizes.")

    matrix1 = np.ascontiguousarray(rm)
    matrix2 = np.ascontiguousarray(tm)
    output = np.ndarray(rm.shape, dtype=np.float64)
    
    p_matrix1 = matrix1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    p_matrix2 = matrix2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    p_output = output.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    c_shape = matrix1.ctypes.shape_as(ctypes.c_int)
    c_dd = ctypes.c_double(dd)
    c_dta = ctypes.c_double(dta)
    
    gamma_index_f(len(matrix1.shape), c_shape, p_matrix1, p_matrix2, p_output, c_dd, c_dta)
    
    return output