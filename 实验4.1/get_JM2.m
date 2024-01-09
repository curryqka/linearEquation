function [J] = get_JM2(x1, x2, x3)
%GET_JM2 此处显示有关此函数的摘要
% 计算Jacobi Matrix
J = [3,x3*sin(x2*x3),x2*sin(x2*x3);
    2*x1, -162*x2 - 81/5, cos(x3);
    -x2*exp(1)^(-x1*x2),-x1*exp(1)^(-x1*x2),20];
end

