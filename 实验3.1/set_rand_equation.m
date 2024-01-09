function [A, b] = set_rand_equation(n)
%SET_RAND_EQUATION 
% 设置自定义矩阵
% Args:dim 维度
% Returns: A,b
A = zeros(n);
x = ones(n, 1);
% for i = 1:n
%     for j = 1:n
%         A(i, j) = 1/(i + j - 1);
%     end
% end
base = [1, 5, 7];
    for i = 1 : n
        if i == 1
            A(i,i:i+1) = base(2:end);
        elseif i == n
                A(i, i-1:i) = base(1:2);
        else
            A(i, i-1:i+1) = base;
        end
    end
b = A * x;
end

