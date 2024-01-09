% --- --- --- --- --- ---定义自由选主元Gauss消去函数 --- --- --- --- --- 
function x = free_gauss(A,b) 
%free_gauss为自由主元高斯消去法
    %% init
    % A dim
    n = size(A,1);
    x = zeros(n,1);
    %l为消元因子
    l = zeros(n,1);
    % e为初等变换矩阵，用于列交换之后将最后的解x的位置更正
    e = eye(n);
    %% 迭代n-1次
    for iter = 1:n-1
        if A(iter, iter) == 0
            disp('系数矩阵奇异,主元为零,完全主元Gauss消元法无法进行');

            break;
        end
        %寻找最大主元
        
        %按模最小或尽可能小选取主元
        MIN = Inf;
        for p=iter:n
            for q=iter:n
                if (abs(A(p,q))<MIN) && (abs(A(p,q))~=0)
                    MIN = abs(A(p,q));
                    index_r = p;
                    index_c = q;
                else
                    index_r = iter;
                    index_c = iter;
                end
            end
        end
        % index_r = index_r+iter-1;
        % index_c = index_c+iter-1;
        % MAX1 = max(abs(A(iter:n,iter:n)));
        % MAX = max(MAX1);
        % [index1,index2] = find(abs(A(iter:n,iter:n))==MAX,1);
        % index1 = index1+iter-1;
        % index2 = index2+iter-1;

        %交换行、列
        if iter ~= index_r
            A([iter,index_r],:) = A([index_r,iter],:);
            b([iter,index_r],:) = b([index_r,iter],:);
        end
        if iter ~= index_c
            A(:,[iter,index_c]) = A(:,[index_c,iter]);
            e(:,[iter,index_c]) = e(:,[index_c,iter]);
        end
        for id=iter+1:n
            l(id) = A(id,iter)/A(iter,iter);
            b(id) = b(id) - l(id)*b(iter);
            A(id, iter:end) = A(id, iter:end) - l(id)*A(iter, iter:end);
            % for j =iter:n
            %     A(id,j) = A(id,j) - l(id)*A(iter,j); 
            % end
        end
    end
    if A(n, n) == 0
        return;
    else
        x(n) = b(n)/A(n,n);
        for iter=n-1:-1:1
            w=0;
            for j = iter+1:n
                w = w + A(iter,j)*x(j);
            end
            x(iter) = (b(iter)-w)/A(iter,iter);
        end
        x=e*x;
    end
end