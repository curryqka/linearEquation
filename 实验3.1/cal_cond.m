function [X,C_1, C_2, C_inf] = cal_cond(A,b)
%CALL_COND 
% 计算矩阵的条件数
%求条件数及理论解
% Args: A,b 系数矩阵与右端项
% Returns: X 理论解; C_1 1-条件数; C_2 2-条件数; C_inf inf-条件数
format long e
disp('线性方程组的精确解:');
X=(A\b); 
disp(X);
fprintf('矩阵A的1-条件数: %f \n',cond(A,1));
fprintf('矩阵A的2-条件数: %f \n',cond(A));
fprintf('矩阵A的无穷-条件数: %f \n',cond(A,inf));
C_1 = cond(A, 1);
C_2 = cond(A);
C_inf = cond(A, inf);

end

