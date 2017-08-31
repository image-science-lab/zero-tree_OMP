y = rand(4,1);
D = rand(4, 16);
l0_resolve = 4;
k0 = 3;
x1 = [1 1 0 0]';

x0 = comp_mex(D, y, x1, k0, l0_resolve);
