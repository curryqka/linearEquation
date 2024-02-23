classdef itermodel
    %ITERMODEL 迭代法的�?
    %   面向解决Hilbert矩阵的类
    
    properties
        H;
        x;
        b;
        n;
        X;
        D;
        L;
        U;
        % w;
        iter_num;
        ee;
    end
    
    methods
        function obj = itermodel(n)
            %ITERMODEL 构�?�此类的实例
            % 构�?�初始化Hilbert矩阵
            % Args: n dimension of the matrix
            obj.H = zeros(n);
            obj.x = ones(n,1);
            obj.b = ones(n,1);
            obj.n = n;
            obj.X = ones(n,1); % 精确�?
            for i = 1:obj.n
                for j = 1:obj.n
                    obj.H(i, j) = 1/(i + j - 1);
                end
            end
            % 构�?�对角矩�?,下三角矩阵和上三角矩�?
            obj.D = diag(diag((obj.H)));
            obj.L = -tril(obj.H, -1);
            obj.U = -triu(obj.H, 1);
            obj.iter_num = 10000;
            % obj.w = w;
            obj.ee = 1e-4;
        end

        % get the parameter b
        function obj = set_matrix(obj, x_)
            %set_matrix 
            % fill the Hilbert
            % Args: x_ 给定�?
            if nargin > 1
                obj.b = obj.H * x_;
            else
                disp('默认b')
                obj.b = obj.H * obj.x;
                disp(obj.b)
            end
        end
        function [loss_1, loss_2, loss_inf] = cal_loss(obj, x)
            error = obj.X - x;
            loss_1 = norm(error, 1);
            loss_2 = norm(error, 2);
            loss_inf = norm(error, inf);
        end

        
        function [L, U] = get_LU(obj)
            %% init L,U
            L = eye(obj.n);
            U = zeros(obj.n);

            %% first element
            U(1, :) = obj.H(1, :);
            L(:, 1) = obj.H(:, 1)/U(1, 1);
            
            %% the step k
            for iter = 2 : obj.n
                U(iter , iter : end) = obj.H(iter, iter : end) - L(iter, :)*U(:, iter : end);
                L(iter + 1:end, iter) = 1/U(iter, iter)*(obj.H(iter + 1:end, iter) ...
                    - L(iter + 1:end, :) * U(:, iter));
            end

        end

        function [L, U] = get_LU_check(obj)
            %   lu decompose
            %   L:下三角矩�?
            %	U:上三角矩�?
            %	A:输入矩阵
            
            A = obj.H;
            L=eye(obj.n);
            L(:,1)=A(:,1)/A(1,1);%L第一列赋�?
            
            U=zeros(obj.n);
            U(1,:)=A(1,:);%U第一行赋�?
            
            for i=2:obj.n
                for j=2:obj.n
                    if i<=j
                        U(i,j)=A(i,j)-sum(L(i,1:i-1).*U(1:i-1,j)');%递推表达�?(1-6)
                    else
                        if U(j,j)==0
                            L(i,j)=0;
                        else
                            L(i,j)=(A(i,j)-sum(L(i,1:j-1).*U(1:j-1,j)'))/U(j,j);%递推表达�?(1-7)
                        end
                    end
                end
            end
        end
        % 高斯迭代(LU分解)
        function x = gauss_method(obj, L, U)
            %% init
            x = ones(obj.n, 1);
            y = ones(obj.n ,1);
            %% LUx = b, Ly = b, Ux = y
            for iter = 1 : obj.n
                y(iter) = obj.b(iter) - L(iter, 1:iter -1)*y(1:iter - 1);
            end
            for iter = obj.n : -1 : 1
                x(iter) = 1/U(iter, iter) * (y(iter) - ...
                    U(iter, iter + 1:end) * x(iter + 1:end));
            end
        end
        
        % 迭代�?
        function rou = cal_rou(obj, mode, w)
            if nargin > 2 && mode == 'S'
                iterMatrix = (obj.D - w*obj.L) \ (( 1- w)*obj.D ...
                    + w*obj.U);
            else
                if mode == 'J'
                    iterMatrix = obj.D \ (obj.L + obj.U);
                elseif mode == 'G'
                    iterMatrix = (obj.D - obj.L) \ obj.U;
              
                else
                    disp('请输入J和G,S中的�?个！')
                end
            end

            rou = max( abs(eig(iterMatrix)) );
        end

        % J�?
        function [x, loss] = jacobi_method(obj)
            %% init
            x = zeros(obj.n, 1);
            loss = zeros(obj.iter_num, 1);
            %% iter matrix J
            B = obj.D \ (obj.L + obj.U);
            f = obj.D \ obj.b;

            iter = 1;
            while iter <= obj.iter_num
            % for iter = 1: obj.iter_num
                [~, ~, loss_inf] = obj.cal_loss(x);
                loss(iter) = loss_inf;
                % save the previous one
                x_ = x;
                
                x = B*x_ + f;

                if norm(x_ - x, 1) < obj.ee
                    disp('J法达到收敛！');
                    disp('迭代次数�?');
                    disp(iter);
                    break;
                end
                iter = iter + 1;
                
            end
            
        end 

        
        % GS�?
        function [x, loss] = gs_method(obj)
            %% init
            x = zeros(obj.n, 1);
            loss = zeros(obj.iter_num, 1);
            %% iter matrix GS
            B = (obj.D - obj.L) \ obj.U;
            f = (obj.D - obj.L) \ obj.b;
            % for iter = 1:obj.iter_num
            iter = 1;
            while iter <= obj.iter_num
                [~, ~, loss_inf] = obj.cal_loss(x);
                loss(iter) = loss_inf;

                x_ = x;

                x = B*x_ + f;
                
                if norm(x_ - x, 1) < obj.ee
                    disp('GS法达到收敛！');
                    disp('迭代次数�?');
                    disp(iter);
                    break;
                end
                iter = iter + 1;

            end
        end
        
        % SOR�?
        function [x, loss] = sor_method(obj, w)
            %% init 
            x = zeros(obj.n, 1);
            loss = zeros(obj.iter_num, 1);
            %% iter matrix SOR
            iter = 1;
            while iter <= obj.iter_num
            % for iter = 1:obj.iter_num
                [~, ~, loss_inf] = obj.cal_loss(x);
                loss(iter) = loss_inf;

                x_ = x;
                x = (obj.D - w*obj.L) \ ( ( (1-w)*obj.D + w*obj.U)*x_ + ...
                    w*obj.b);

                if norm(x_ - x, 1) < obj.ee
                    disp('SOR法达到收敛！');
                    disp('迭代次数�?');
                    disp(iter);
                    break;
                end
                iter = iter + 1;
            end
        end

        % draw the loss curve
        function get_losscurve(obj, error1, name1,...
                error2, name2, error3, name3)
            figure;
            maxIter = obj.iter_num;
            lw = 2.5;
            if mod(nargin - 1, 2) == 0
                if nargin - 1 == 2
                    % plot(1 : maxIter, error1, LineWidth=lw);
                    % 2018 can not LineWidth=lw;
                    plot(1 : maxIter, error1, 'LineWidth', lw);
                    hold on
                    ax = gca;
                    ax.LineWidth = 1.5;
                    ax.FontName = '微软雅黑';
                    ax.FontSize = 14;
                    grid on
                    set(gca,'GridLineStyle','--')
                    xlabel('迭代次数');
                    ylabel('误差');
                    title('线�?�方程组迭代法误差曲�?');
                    legend(name1);
                elseif nargin - 1 == 4
                    %plot(1 : maxIter, error1, LineWidth=lw);
                    plot(1 : maxIter, error1, 'LineWidth', lw)
                    hold on
                    plot(1 : maxIter, error2, 'LineWidth', lw);
                    ax = gca;
                    ax.LineWidth = 1.5;
                    ax.FontName = '微软雅黑';
                    ax.FontSize = 14;
                    grid on
                    set(gca,'GridLineStyle','--')
                    xlabel('迭代次数');
                    ylabel('误差');
                    title('线�?�方程组迭代法误差曲�?');
                    legend(name1, name2);
                elseif nargin - 1 == 6
                    plot(1 : maxIter, error1, 'LineWidth', lw);
                    hold on
                    plot(1 : maxIter, error2, 'LineWidth', lw);
                    plot(1 : maxIter, error3, 'LineWidth', lw);
                    ax = gca;
                    ax.LineWidth = 1.5;
                    ax.FontName = '微软雅黑';
                    ax.FontSize = 14;
                    grid on
                    set(gca,'GridLineStyle','--')
                    xlabel('迭代次数');
                    ylabel('误差');
                    title('线�?�方程组迭代法误差曲�?');
                    legend(name1, name2, name3);
                else
                    disp('函数输出参数数目�?2,4,6')
                end
            else
                disp('输入函数参数不为偶数!')
            end
           
            % legend('Jacobi迭代�?', 'Gauss-Seidel迭代�?', 'SOR迭代�?');
        
        end
    end
end

