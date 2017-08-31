% Script to train dictionaries with various methods.
clear all;
close all;

addpath(genpath('tools'));

%% Parameters

% (h, w, t) is the size of video patch. Change it here to reflect in your
% training algorithm
h = 8; w = 8; t = 16;
R = get_video_resizer(h/2, h, t/2, t);  % Undersampling matrix.
N1 = 512;                               % Low resolution dictionary.
l0_resolve = 16;                        % Number of high resolution atoms per atom.
k1 = 2;                                 % Low resolution sparsity.
k0 = 4;                                 % High resolution sparsity.
niters = 10;                            % Number of K-SVD iterations.
save_dir = 'data';
ndata = 100000;

W = R;                                  % Downsampler.
U = 8*R';                               % Upsampler.

% Create the directory to save output.
mkdir(save_dir);

%% Data path.
% Replace '/path/to/your/data' with the appropriate path. Remember that the
% training data is a matrix of dimensions (h*w*t)X(number of data points)
fprintf('Loading data\n');
tic; load('/path/to/your/data'); data_time = toc;
fprintf('Loaded data. Time take = %f\n', data_time);
data = Y; clear Y;
return;

% Use fewer data points for training
data = data(:, 1:ndata);

%% Call master script now
train_dictionaries;
