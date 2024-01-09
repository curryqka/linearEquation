clear;
clc;
format long;

%方法选择
n=input('矩阵A阶数:n=');
disp('选取求解方式');
disp('1 顺序Gauss消元法,2 列主元Gauss消元法,3 完全选主元Gauss消元法,4 模最小或近可能小的元素作为主元');
a=input('求解方式序号:');

%赋值A和b
A=zeros(n,n);
b=zeros(n,1);
for i=1:n    
   A(i,i)=6;
     if i>1
          A(i,i-1)=8;
     end
     if i<n
          A(i,i+1)=1;
     end
end
for i=1:n   
    for j=1:n
    b(i)=b(i)+A(i,j);
    end
end
disp('给定系数矩阵为:A');
disp(A);
disp('右端向量为:');
disp(b);

%求条件数及理论解
disp('线性方程组的精确解:');
X=(A\b); 
fprintf('矩阵A的1-条件数: %f \n',cond(A,1));
fprintf('矩阵A的2-条件数: %f \n',cond(A));
fprintf('矩阵A的无穷-条件数: %f \n',cond(A,inf));
    
 

%顺序Gauss消元法  
if a==1  
    A1=A;b1=b;
    for k=1:n
        if A1(k,k)==0
            disp('主元为零,顺序Gauss消元法无法进行');
            break
        end
    fprintf('第%d次消元所选取的主元:%g\n',k,A1(k,k))
    disp('此次消元后系数矩阵为:');
    disp(A1);
        for p=k+1:n
            l=A1(p,k)/A1(k,k);
            A1(p,k:n)=A1(p,k:n)-l*A1(k,k:n);
            b1(p)=b1(p)-l*b1(k);
        end
    end
x1(n)=b1(n)/A1(n,n);
for k=n-1:-1:1
   for w=k+1:n        
        b1(k)=b1(k)-A1(k,w)*x1(w);
   end
    x1(k)=b1(k)/A1(k,k);
end
    disp('顺序Gauss消元法解为:');
    disp(x1);
    disp('所求解与精确解之差的无穷-范数为');
    norm(x1-X,inf) 
end

%列主元Gauss消元法
if a==2
    A2=A;b2=b;
    for k=1:n
        [max_i,max_j]=find(A2(:,k)==max(abs(A2(k:n,k))));
        if max_i~=k
            A2_change=A2(k,:);
            A2(k,:)=A2(max_i,:);
            A2(max_i,:)=A2_change;
            b2_change=b2(k);
            b2(k)=b2(max_i);
            b2(max_i)=b2_change;
        end
        if A2(k,k)==0
            disp('主元为零,列主元Gauss消元法无法进行');
            break
        end
    fprintf('第%d次消元所选取的主元:%g\n',k,A2(k,k))
    disp('此次消元后系数矩阵为:');
    disp(A2);
        for p=k+1:n
            l=A2(p,k)/A2(k,k);
            A2(p,k:n)=A2(p,k:n)-l*A2(k,k:n);
            b2(p)=b2(p)-l*b2(k);
        end
    end
x2(n)=b2(n)/A2(n,n);
for k=n-1:-1:1
   for w=k+1:n        
        b2(k)=b2(k)-A2(k,w)*x2(w);
   end
    x2(k)=b2(k)/A2(k,k);
end
    disp('列主元Gauss消元法解为:');
    disp(x2);
    disp('所求解与精确解之差的无穷-范数为');
    norm(x2-X,inf) 
end

%完全选主元Gauss消元法
if a==3  
    A3=A;b3=b;
    for k=1:n
        VV=eye(n);
        [max_i,max_j]=find(A3(k:n,k:n)==max(max(abs(A3(k:n,k:n)))));
        if numel(max_i)==0
            [max_i,max_j]=find(A3(k:n,k:n)==-max(max(abs(A3(k:n,k:n)))));
        end
        W=eye(n);
        W(max_i(1)+k-1,max_i(1)+k-1)=0;
        W(k,k)=0;
        W(max_i(1)+k-1,k)=1;
        W(k,max_i(1)+k-1)=1;
        V=eye(n);
        V(k,k)=0;
        V(max_j(1)+k-1,max_j(1)+k-1)=0;
        V(k,max_j(1)+k-1)=1;
        V(max_j(1)+k-1,k)=1;
        A3=W*A3*V;
        b3=W*b3;
        VV=VV*V;

        if A3(k,k)==0
            disp('主元为零,完全选主元Gauss消元法无法进行');
            break
        end
    fprintf('第%d次消元所选取的主元:%g\n',k,A3(k,k))
    disp('此次消元后系数矩阵为:');
    disp(A3)
        for p=k+1:n
            l=A3(p,k)/A3(k,k);
            A3(p,k:n)=A3(p,k:n)-l*A3(k,k:n);
            b3(p)=b3(p)-l*b3(k);
        end
    end
x3(n)=b3(n)/A3(n,n);
for k=n-1:-1:1
   for w=k+1:n        
        b3(k)=b3(k)-A3(k,w)*x3(w);
   end
    x3(k)=b3(k)/A3(k,k);
end
    disp('完全选主元Gauss消元法解为:');
    disp(x3);
    disp('所求解与精确解之差的无穷-范数为');
    norm(x3-X,inf) 
end

%模最小或近可能小的元素作为主元
if a==4  
    A4=A;b4=b;
    for k=1:n
        AA=A4;
        AA(AA==0)=NaN;
        [min_i,j]=find(AA(k:n,k)==min(abs(AA(k:n,k))));
        if numel(min_i)==0
            [min_i,j]=find(AA(k:n,k)==-min(abs(AA(k:n,k:n))));
        end
        W=eye(n);
        W(min_i(1)+k-1,min_i(1)+k-1)=0;
        W(k,k)=0;
        W(min_i(1)+k-1,k)=1;
        W(k,min_i(1)+k-1)=1;
        A4=W*A4;
        b4=W*b4;
        if A4(k,k)==0
            disp('主元为零,模最小Gauss消元法无法进行');
            break
        end
    fprintf('第%d次消元所选取的主元:%g\n',k,A4(k,k))
    %A4
        for p=k+1:n
            l=A4(p,k)/A4(k,k);
            A4(p,k:n)=A4(p,k:n)-l*A4(k,k:n);
            b4(p)=b4(p)-l*b4(k);
        end
    end
x4(n)=b4(n)/A4(n,n);
for k=n-1:-1:1
   for w=k+1:n        
        b4(k)=b4(k)-A4(k,w)*x4(w);
   end
    x4(k)=b4(k)/A4(k,k);
end
    disp('模最小Gauss消元法解为:');
    disp(x4);
    disp('所求解与精确解之差的无穷-范数为');
    norm(x4-X,inf) 
end