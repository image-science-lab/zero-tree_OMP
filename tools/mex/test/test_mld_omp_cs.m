addpath(genpath('../'));
addpath(genpath('../../'));

D = rand(128, 1024);
norms =diag(sqrt(D'*D));
D = bsxfun(@rdivide, D, norms');
y = rand(128,1);
k = 5;
thres = 1e-4;
blksize = 64;

A = sprand(64, 128, 0.5);

tic; x_c = mld_omp_cs_c(D, A*y, A, blksize, thres); toc

Dmod = A*D;
norms = diag(sqrt(Dmod'*Dmod));
Dmod = bsxfun(@rdivide, Dmod, norms');

tic;
x_b = mld_omp(Dmod, A*y, blksize, thres)./(norms*ones(1, size(y, 2)));
toc

