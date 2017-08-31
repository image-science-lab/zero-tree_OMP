%% Script to test if sdot is fast.

clear all;
close all;

m = 512;
n = 8192;
k = 64;

D = randn(m, n);
y = randn(m, 1);
sup = randi(m, 1, k);

tic; Dsub = D(:, sup); x1 = Dsub'*y; time1 = toc
tic; x2 = sdot(D, y, uint32(sup)); time2 = toc
