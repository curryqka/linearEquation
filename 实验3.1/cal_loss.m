function [loss_1, loss_2, loss_inf] = cal_loss(x_, x)
%CAL_LOSS 
% calculate the error
format long e
error = x_ - x;
loss_1 = norm(error, 1);
loss_2 = norm(error, 2);
loss_inf = norm(error,inf);
end

