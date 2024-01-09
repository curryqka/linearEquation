function [x] = col_gauss(A,b)
%COL_GAUSS 
% 列主元高斯消去
% 消去
    for iter = 1 : length(b)-1
        % find the max aij (column)
        if A(iter, iter) == 0
            disp('系数矩阵奇异,主元为零,列主元Gauss消元法无法进行');

            break;
        end
        [id,~] = find( A(:,iter) == max(abs(A(iter:end, iter))) );
        if id ~= iter
        % id = id + iter - 1;
        % exchange
        % tmp = A(iter,iter:end);
        % A(iter,iter:end) = A(id,iter:end);
        % A(id,iter:end) = tmp;
            tmp = A(iter,:);
            A(iter,:) = A(id,:);
            A(id,:) = tmp;

            tmp = b(iter);
            b(iter) = b(id);
            b(id) = tmp;
        end
        
        for row = iter + 1 : length(b)
            % 消元
            L = -A(row, iter) / A(iter, iter);
            % 相减
            A(row, iter:end) = A(row, iter:end) + L*A(iter, iter:end);
            b(row) = b(row) + L*b(iter);
        end
    end
    % 回代
    x = zeros(length(b),1);
    for iter = length(b):-1:1
        if A(length(b), length(b)) == 0
            break;
        end
        if iter == length(b) % 第n个分量
            x(iter) = b(iter)/A(iter, iter);
        else
            x(iter) = (b(iter) - A(iter, iter + 1 : end) * x(iter + 1 : end) ) / A(iter, iter);
        end
    end
end

