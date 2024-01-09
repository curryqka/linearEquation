function [save_X] = newton_2(init_value, delta, iter_num)
%NEWTON_2 此处显示有关此函数的摘要
% solve the equation(2)
save_X = [];

x1 = init_value(1); x2 = init_value(2); x3 = init_value(3);
init_X = init_value;
iter = 0;
%% 牛顿迭代法
tic;
while iter <= iter_num
    % set equation
    iter = iter + 1;

    F = get_F2(x1, x2, x3);
    % jacobi_M=[12,-2*x2,-4;2*x1,10,-1;0, 3*x2^2,10];
    jacobi_M = get_JM2(x1, x2, x3);
    % dx=-jacobi_M\F; %matlab中默认是列主元消去解方程组     
    % x=x+dx;
    % error=max(abs(dx));

    delta_X = jacobi_M \ (-F);
    X = init_X + delta_X;
    error = norm(X - init_X,1);
    init_X = X;

    % x1 = X(1); x2 = X(2); x3 = X(3);
    [x1, x2, x3] = get_new_X(X);

    if error < delta
        
        disp('Newtons法满足设定误差!');
        disp('The IterNum is:');
        disp(iter)
        break;
    end
end
toc;

save_X = [save_X, X];
%% 拟牛顿迭代法
% init value
x1 = init_value(1); x2 = init_value(2); x3 = init_value(3);
init_X = [x1; x2; x3];
% init_A = F'(X0)
% init_A = [12, -2*init_X(2), -4;
%     2*init_X(1), 10, -1;
%     0, 3*init_X(2)^2, 10];

init_A = get_JM2(x1, x2, x3);
iter = 0;
tic
while iter <= iter_num
    iter = iter + 1;
    F = get_F2(x1, x2, x3);
    
    delta_X = -init_A \ F;
    
    X = init_X + delta_X;
    
    error = norm(X - init_X, 1);

    init_X = X;
    [x1, x2, x3] = get_new_X(X);
    

    if error < delta
        
        disp('拟Newtons法满足设定误差!');
        disp('The IterNum is:');
        disp(iter)
        break;
    end
    % upgrate F
    F_1 = get_F2(x1, x2, x3);
    Y = F_1 - F;
    A = init_A + (Y - init_A * delta_X)*delta_X'/(delta_X'*delta_X);
    
    init_A = A;
    
end
toc
save_X = [save_X, X];
end

