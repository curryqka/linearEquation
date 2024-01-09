%% clc gauss31_main.m
clc;clear;close all
format long e;
%% 系数矩阵f 右端项b
% test: 
num_matrix = 50;
% [A,b] = set_equation(num_matrix);
[A, b] = set_rand_equation(num_matrix);
[X, C1, C2, Cinf] = cal_cond(A, b)
X = ones(num_matrix,1);

%% 顺序gauss
x = order_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('顺序高斯x:');
disp(x);

%% 列主元gauss
x = col_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('列高斯x:');
disp(x);

%% 完全主元gauss
x = all_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('完全高斯x:');
disp(x);

%% 自由主元gauss (e.g. min)
x = free_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('自由高斯x:');
disp(x)