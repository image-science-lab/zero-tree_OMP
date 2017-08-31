% Script for recovering video using Hitomi's method.

clear all;
close all;

addpath(genpath('../utils'));

% Load data
load('data/two_scale_dict.mat');
load('data/single_scale_dict.mat');

load('data/test_sharpner.mat');

[H, W, T] = size(video);
h = 8; w = 8; t = 16;
sh = 8; sw = 8; st = 8;
R = sparse(R);

% Construct the sampling operator.
code_block = get_code_block(h, w, t/2);
code_block = repmat(code_block, 1, 1, 2);

nz_idx = find(code_block == 1);
nlen = length(nz_idx);

A = sparse(1:nlen, nz_idx, ones(nlen, 1), nlen, h*w*t);

% Get video patches.
Y = vid2mat(video, h, w, t, sh, sw, st);

% Recover using OMP.
fprintf('Recovering using OMP\n');
tic; X = omp_cs_c(D, A*Y, A, k1, 1e-4); omp_time = toc;
Y_omp = D*X;
snr_omp = rsnr(Y, Y_omp);
vid_omp = mat2vid(Y_omp, h, w, t, sh, sw, st, H, W, T);

% Recover using zero tree OMP.
fprintf('Recovering using zero tree OMP\n');
tic; [~, X0] = comp_cs_c(D1, D0, A*Y, 8*R', A, k1, k0, 1e-4);
comp_time = toc;
Y_comp = D0*X0;
snr_comp = rsnr(Y, Y_comp);
vid_comp = mat2vid(Y_comp, h, w, t, sh, sw, st, H, W, T);

% Save videos.
create_video(uint8(vid_omp), 'omp_cs', 30, false);
create_video(uint8(vid_comp), 'comp_cs', 30, false);

fprintf('OMP SNR: %.2f; zero tree OMP SNR: %.2f\n', snr_omp, snr_comp);
fprintf('OMP time: %.2f; zero tree OMP time: %.2f\n', omp_time, comp_time);
fprintf('Speed up: %.2f\n', omp_time/comp_time);
