addpath('../../../scale_ksvd');
addpath('../../');
addpath(genpath('../../../../external_code/ompbox'));
A = rand(128, 1024);
norms =diag(sqrt(A'*A));
A = bsxfun(@rdivide, A, norms');
y = rand(128,1);
k = 5;
thres = 1e-4;

tic; x_c = omp_c(A, y, k, thres); toc
tic;
x_m = zeros(size(A, 2), size(y, 2));
    for idx = 1:1:size(y,2)
        x_m(:, idx) = my_omp(A, y(:,idx), k, 1e-4);
    end

toc

tic;
x_b = solveX(y, A, k);
toc

