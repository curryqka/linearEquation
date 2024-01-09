function [jacobi_M] = get_JM1(x1, x2, x3)
%GET_JM 此处显示有关此函数的摘要
% 计算Jacobi矩阵
jacobi_M=[12,-2*x2,-4;2*x1,10,-1;0 3*x2^2,10];
end

