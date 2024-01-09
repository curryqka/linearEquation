function F = get_F1(x1, x2, x3)
% upgrate the value
    f1=12*x1-x2^2-4*x3-7;     
    f2=x1^2+10*x2-x3-11;     
    f3=x2^3+10*x3-8;
    % x1 = x1;
    % x2 = x2;
    % x3 = x3;
    % f1 = subs(F1);
    % f2 = subs(F2);
    % f3 = subs(F3);
    
    F=[f1,f2,f3]';
end