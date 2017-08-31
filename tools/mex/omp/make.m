% Make file for testing subroutines. 
output = '../bin';
mkdir(output);
%libpath = '/home/vishwanath/lib';
libpath = '../../lib';

sources = {'../utils.c', '../linalg.c', 'omp_core.c', '../sparse.c'};
args = {'-llapack', '-lblas', '-lgfortran', sprintf('-L%s', libpath), '-I../', '-I.', '-largeArrayDims'};

files = {'dxmul.c', 'comp_c3.c', 'comp_cs_c.c', 'omp_c.c', 'omp_cs_c.c', 'comp_c.c', 'mld_omp_c.c', 'mld_omp_cs_c.c'};
%files = {'comp_cs_c.c'};

for idx = 1:length(files)
	fprintf('Compiling %s\n', files{idx});
	mex(files{idx}, sources{:}, args{:});
end
