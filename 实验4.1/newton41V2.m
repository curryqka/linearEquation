%% clean the room
clc;clear;close 

%% init x
% test example: given init_X, use function : get_experiment_ans
init_X = [2;2;2];
[save_X_1, save_X_2, save_iter] = get_experiment_ans(init_X)

% init_X = [0;0;0];
% [save_X_1, save_X_2, save_iter] = get_experiment_ans(init_X)

function [save_X_1, save_X_2, save_iter] = get_experiment_ans(init_X)
%% property
delta = 0.5 * 1e-5;
iter_num = 1000;

% init value
% x1 = 1; x2 = 1; x3 = 1;
% init_X = [1;1;1];

% 存储解的矩阵:[3, 2] column 1: newton method; column 2: broyden newton method
save_X_1 = zeros(3,2);
save_X_2 = zeros(3,2);
% 存储迭代次数的矩阵:[2, 2] row: equation, column: method
% e.g. save_iter(1,1) 第一个方程,newton法的迭代次数(第一种方法)
save_iter = ones(2);
%% Solve
% 类实体 newton : solve equation(1)
equation1 = NewtonMethod(myfun1(), init_X, delta, iter_num);

% 求解类实体:newton1 newton method
[newton1, iter_newton1] = equation1.newton();
save_X_1(:, 1) = newton1.X;

% 求解类实体:b_newton1 broyden newton method
[b_newton1,iter_b_newton1] = equation1.b_newtons();
save_X_1(:, 2) = b_newton1.X;

save_iter(1,:) = [iter_newton1, iter_b_newton1];
% 类实体 newton : solve equation(2)
equation2 = NewtonMethod(myfun2(), init_X, delta, iter_num);

% 求解类实体:newton2 newton method
[newton2, iter_newton2] = equation2.newton();
save_X_2(:, 1) = newton2.X;

% 求解类实体:b_newton2 broyden newton method
[b_newton2, iter_b_newton2] = equation2.b_newtons();
save_X_2(:, 2) = b_newton2.X;

save_iter(2,:) = [iter_newton2, iter_b_newton2];
end

function f = myfun1()

    syms x1 x2 x3
    f1=12*x1-x2^2-4*x3-7;     
    f2=x1^2+10*x2-x3-11;     
    f3=x2^3+10*x3-8;
    f = [f1;f2;f3];
end

function f = myfun2()

    syms x1 x2 x3
    f1=3*x1-cos(x2*x3)-0.5;     
    f2=x1^2-81*(x2+0.1)^2+sin(x3)+1.06;     
    f3=exp(1)^(-x1*x2)+20*x3+(10*pi-3)/3;
    f=[f1,f2,f3]';
end


% syms x1 x2 x3
% f1=12*x1-x2^2-4*x3-7;     
% f2=x1^2+10*x2-x3-11;     
% f3=x2^3+10*x3-8;
% 
% x1 = 1; x2 = 1; x3 = 1;
% F = get_F1(f1, f2, f3, x1, x2, x3)

% test jacobian
% syms x1 x2
% J = jacobian(myfun(x1, x2), [x1 x2])
% 
% 
% JM = get_JM(J, 2, 2)
% 
% function JM = get_JM(J,x1,x2)
% 
% x1 = x1;
% x2 = x2;
% JM = subs(J);
% end
% function y = myfun(x1, x2)
% y(1) = x1^2 + x2^2;
% y(2) = x1 + x2;
% end

