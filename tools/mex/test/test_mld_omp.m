addpath('../../');
addpath('../omp');

cd ../omp; make; cd ../test/;

D = randn(64, 1024);
norms = diag(sqrt(D'*D));
D = bsxfun(@rdivide, D, norms');

k = 5;
blksize = 64;
thres = 1e-4;

y = randn(64, 2);

tic; x1 = mld_omp(D, y, blksize, thres); toc
tic; x2 = mld_omp_c(D, y, blksize, thres); toc

