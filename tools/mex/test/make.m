% Make file for testing subroutines. 

sources = {'../utils.c', '../linalg.c', '../sparse.c'};

args = {'-llapack', '-lblas', '-lgfortran', '-L../../lib', '-I../', '-largeArrayDims'};

mex('test_spdotm.c', args{:}, sources{:});
