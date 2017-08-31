% Make file for testing subroutines. 
output = '../bin';
mkdir(output);

% Create blas and lapack folders
lapacklib = fullfile(matlabroot, 'extern', computer('arch'), 'microsoft', 'libmwlapack.lib');
blaslib = fullfile(matlabroot, 'extern', computer('arch'), 'microsoft', 'libmwblas.lib');

sources = {'../utils.c', '../linalg.c', 'omp_core.c', '../sparse.c'};
args = {'-I../', '-I.', lapacklib, blaslib};

files = {'comp_c3.c', 'comp_cs_c.c', 'omp_c.c', 'omp_cs_c.c', 'comp_c.c', 'mld_omp_c.c', 'mld_omp_cs_c.c'};

for idx = 1:length(files)
	fprintf('Compiling %s\n', files{idx});
	mex('-largeArrayDims', files{idx}, sources{:}, args{:});
end
