%% clean the room
clc;clear;close 

%% property
delta = 0.5 * 1e-5;
iter_num = 1000;

% init value
x1 = 1; x2 = 1; x3 = 1;
init_X = [1;1;1];

result = newton_1(init_X, delta, iter_num)

result = newton_2(init_X, delta, iter_num)


newton = NewtonMethod(myfun(), [1;1;1], delta, iter_num)

newton1 = newton.newton()

newton1.X

newton2 = newton.b_newtons()

newton2.X

function f = myfun()

    syms x1 x2 x3
    f1=12*x1-x2^2-4*x3-7;     
    f2=x1^2+10*x2-x3-11;     
    f3=x2^3+10*x3-8;
    f = [f1;f2;f3];
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

