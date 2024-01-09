function [A,b] = set_equation(n)
    %SET_MATRIX 
    % args: n 矩阵维度
    % returns: A 生成矩阵 b
    % [6 1,...,
    %    8,6,1,...
    %       
    %       8,6,1
    %         8,6]
    % [7, 15, 15,...,14]
    A = zeros(n);
    base = [8, 6, 1];
    for i = 1 : n
        if i == 1
            A(i,i:i+1) = base(2:end);
        elseif i == n
                A(i, i-1:i) = base(1:2);
        else
            A(i, i-1:i+1) = base;
        end
    end

    b = zeros(n, 1);
    b = sum(A,2);
        
end

