classdef NewtonMethod
    %NEWTONMETHOD 此处显示有关此类的摘要
    % 牛顿法解方程
    
    properties
        % def init value
        init_X;
        % def function
        F;
        % def delta : accuracy
        delta;
        % def number of iteration
        iter_num;
        % def obj solution
        X;
        % def init value
        init_value;
        % def sym
        symbolX;
    end
    
    methods
        function obj = NewtonMethod(F, init_X, delta, iter_num)
            %NEWTONMETHOD 构造此类的实例
            %   此处显示详细说明
            obj.F = F;
            obj.init_X = init_X;
            obj.delta = delta;
            obj.iter_num = iter_num;
            obj.X = init_X;
            % save the init value
            obj.init_value = init_X;
            obj = obj.flush_state();
        end
        
        function obj = flush_state(obj)
            syms x1 x2 x3
            obj.symbolX = [x1;x2;x3];
            obj.init_X = double(obj.init_value);
            obj.X = double(obj.init_value);
        end
        function result = get_new_F(obj)
            % get_F 此处显示有关此方法的摘要
            %   此处显示详细说明:计算F(X)
            % x1 = x1;
            % x2 = x2;
            % x3 = x3;
            x1 = obj.symbolX(1);
            x2 = obj.symbolX(2);
            x3 = obj.symbolX(3);
            [t1, t2, t3] = obj.get_new_X();
            f1 = subs(obj.F(1), {x1, x2, x3}, {t1, t2, t3});
            f2 = subs(obj.F(2), {x1, x2, x3}, {t1, t2, t3});
            f3 = subs(obj.F(3), {x1, x2, x3}, {t1, t2, t3});
            result = [f1;f2;f3];
            result = double(result);
        end

        function JM = get_JM(obj)
            % get_F 此处显示有关此方法的摘要
            %   此处显示详细说明:计算F'(X)
            x1 = obj.symbolX(1);
            x2 = obj.symbolX(2);
            x3 = obj.symbolX(3);
            J = jacobian(obj.F, [x1, x2, x3]);
            [t1, t2, t3] = obj.get_new_X();
            JM = subs(J,{x1,x2,x3},{t1,t2,t3});
            JM = double(JM);
        end

        function [x1, x2, x3] = get_new_X(obj)
            x1 = obj.X(1);
            x2 = obj.X(2);
            x3 = obj.X(3);
        end

        function [obj, iter] = newton(obj)
            % newton method to solve equation
            obj = obj.flush_state();
            
            iter = 0;
            tic;
            while iter <= obj.iter_num
                iter = iter + 1;
                if mod(iter, 5) == 0
                    disp('newton---Processing Iter');
                    disp(iter);
                end
                % 计算F(x)和F'(x)
                ff = obj.get_new_F();
                jacobi_M = obj.get_JM();
                
                delta_X = jacobi_M \ (-ff);
                
                obj.X = obj.init_X + delta_X;
                % delta_X
                if norm(delta_X, 1) < obj.delta
                    disp('Newton Method 满足预设误差!')
                    disp('The IterNum is:');
                    disp(iter);
                    break;
                end
                obj.init_X = obj.X;

            end
            toc;   
        end

        function [obj, iter] = b_newtons(obj)
            % broyden newton method to solve equation
            iter = 0;
            obj = obj.flush_state();
            % obj.X
            init_A = obj.get_JM();
            tic;
            while iter <= obj.iter_num
                iter = iter + 1;
                if mod(iter,5) == 0
                    disp('Broyden newton---Processing Iter：');
                    disp(iter);
                end
                ff = obj.get_new_F();
                
                % 计算Ak,计算F(X)
                delta_X = -init_A \ ff;
                
                obj.X = obj.init_X + delta_X;

                if norm(delta_X, 1) < obj.delta
                    disp('Broyden Newton Method 满足预设误差!')
                    disp('The IterNum is:');
                    disp(iter);
                    break;
                end
                obj.init_X = obj.X;
                % upgrate F
                F_1 = obj.get_new_F();
                Y = F_1 - ff;
                A = init_A + (Y - init_A * delta_X)*delta_X'/(delta_X'*delta_X);
                
                init_A = A;
            end
            toc;
            
        end
    end
end

