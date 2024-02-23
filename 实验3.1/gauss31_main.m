%% clc gauss31_main.m
clc;clear;close all
format long e;
%% 系数矩阵f 右端项b
% test: 
n=input('矩阵A阶数:n=');
num_matrix = n;
[A, b] = set_equation(num_matrix);
% [A, b] = set_rand_equation(num_matrix);
[X, C1, C2, Cinf] = cal_cond(A, b)
X = ones(num_matrix,1);

%方法选择

disp('选取求解方式');
disp('1 顺序Gauss消元法,2 列主元Gauss消元法,3 完全选主元Gauss消元法,4 模最小或近可能小的元素作为主元');
a=input('求解方式序号:');

%% 顺序gauss
if a == 1
x = order_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('顺序高斯x:');
disp(x);

elseif a == 2
%% 列主元gauss
x = col_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('列高斯x:');
disp(x);

elseif a == 3
%% 完全主元gauss
x = all_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('完全高斯x:');
disp(x);

elseif a == 4
%% 自由主元gauss (e.g. min)
x = free_gauss(A,b);
[~, ~, loss_inf] = cal_loss(X, x)
disp('自由高斯x:');
disp(x)

else
    disp('请输入1,2,3,4中的一个!')
end