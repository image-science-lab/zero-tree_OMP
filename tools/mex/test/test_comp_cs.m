addpath(genpath('../'));
addpath(genpath('../../'));

D0 = rand(16, 64);
norms =diag(sqrt(D0'*D0));
D0 = bsxfun(@rdivide, D0, norms');

D1 = rand(4, 16);
norms =diag(sqrt(D1'*D1));
D1 = bsxfun(@rdivide, D1, norms');

y = rand(16,1);
k1 = 3;
k0 = 3;
thres = 1e-4;
hf = 1;

A = sprand(10, 16, 0.2);
U = sprand(16, 4, 0.2);

tic; [x_c1, x_c0] = comp_cs_c(D1, D0, A*y, U, A, k1, k0, thres, hf); toc
x_c0

tic; [x_m1, x_m0] = comp_cs(D1, D0, A*y, U, A, k1, k0, thres, hf); toc
x_m0

