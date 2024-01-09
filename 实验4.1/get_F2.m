function [F] = get_F2(x1, x2, x3)
%GET_F2 此处显示有关此函数的摘要
% get the equation(2) value
f1=3*x1-cos(x2*x3)-0.5;     
f2=x1^2-81*(x2+0.1)^2+sin(x3)+1.06;     
f3=exp(1)^(-x1*x2)+20*x3+(10*pi-3)/3;
F=[f1,f2,f3]';
end

