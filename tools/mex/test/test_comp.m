addpath('../../../scale_ksvd');
addpath('../../');
addpath(genpath('../../../../external_code/ompbox'));
addpath(genpath('../'));
A1 = rand(16, 64);
norms = diag(sqrt(A1'*A1));
A1 = bsxfun(@rdivide, A1, norms');

A0 = rand(64, 1024);
norms = diag(sqrt(A0'*A0));
A0 = bsxfun(@rdivide, A0, norms');

y = rand(64,1);
k1 = 5; k0 = 5;
thres = 1e-4;

R = sparse(get_resizer(4, 8, false));

tic; [x_c1, x_c2] = comp_c(A1, A0, y, R, k1, k0, thres); t_comp = toc;

tic;
x_m1 = solveX(R*y, A1, k1);
x_m0 = solveXConstrained(y, A0, x_m1, k0, 16);
toc

tic; x_c = omp_c(A0, y, k0, thres); t_omp = toc;
tic; x_m = solveX(y, A0, k0); toc

tau = t_omp/t_comp
