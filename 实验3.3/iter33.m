clc;clear;close all


% test example : given dim, use class itermodel
dim = 6;
% dim = 8;
iterModel = itermodel(dim);
%% matlab 要返回类才会更新
iterModel = iterModel.set_matrix();
% disp('精确解：')
% X = inv(iterModel.H) * iterModel.b;

%% 定义需要存储的数据
save_loss_inf = [];
save_x = [];
%% LU分解
[L, U] = iterModel.get_LU();
[L, U] = iterModel.get_LU_check();
disp('LU分解：')
x = iterModel.gauss_method(L, U)
[~, ~, loss_inf] = iterModel.cal_loss(x)

save_loss_inf = [save_loss_inf; loss_inf];
save_x = [save_x, x];

%% J法
rou_J = iterModel.cal_rou('J')
disp('jacobi：')
[x, error_jacobi] = iterModel.jacobi_method();
[~, ~, loss_inf] = iterModel.cal_loss(x)
% save_loss_inf = [save_loss_inf; loss_inf];

save_x = [save_x, x];

%% GS法
rou_GS = iterModel.cal_rou('G')
disp('高斯塞达尔：')
[x, error_gauss_seidel] = iterModel.gs_method();
[~, ~, loss_inf] = iterModel.cal_loss(x)

save_loss_inf = [save_loss_inf; loss_inf];
save_x = [save_x, x];
%% SOR法
rou_SOR = iterModel.cal_rou('S', 0.5)
disp('sor：')
[x, error_sor] = iterModel.sor_method(0.5);
[~, ~, loss_inf] = iterModel.cal_loss(x)
% draw the loss fig
if max(error_jacobi) > 100
    iterModel.get_losscurve(error_gauss_seidel, 'GS',error_sor, 'SOR')
else
    iterModel.get_losscurve(error_jacobi, 'J', error_gauss_seidel, 'GS',error_sor,'SOR')
end

save_loss_inf = [save_loss_inf; loss_inf];
save_x = [save_x, x];
%% draw loss_inf
figure
lw = 1.5;
t = 1:length(save_loss_inf);
bar(t, save_loss_inf, LineWidth=lw);
hold on

for i = 1:length(save_loss_inf)
    
    text(t(i),save_loss_inf(i),num2str(save_loss_inf(i),'%g%'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
 
end


ax = gca;
ax.LineWidth = 1.5;
ax.FontName = '微软雅黑';
ax.FontSize = 14;
grid on
set(gca,'GridLineStyle','--')

set(gca, 'XTick', 1:3, 'XTickLabel',{'LU分解','GS法','SOR法'})
xlabel('不同的方法');
ylabel('误差');
% title('不同w的迭代法误差曲线');
legend(['维度',num2str(dim)]);


%% SOR comparsion
save_error_sor = [];
% w = 0.5
rou_SOR = iterModel.cal_rou('S', 0.5)
disp('sor w=0.5：')
[x, error_sor] = iterModel.sor_method(0.5);
[~, ~, loss_inf] = iterModel.cal_loss(x)
save_error_sor = [save_error_sor, error_sor];

% w = 0.8
rou_SOR = iterModel.cal_rou('S', 0.8)
disp('sor w=0.8：')
[x, error_sor] = iterModel.sor_method(0.8);
[~, ~, loss_inf] = iterModel.cal_loss(x)
save_error_sor = [save_error_sor, error_sor];

% w = 1.3
rou_SOR = iterModel.cal_rou('S', 1.3)
disp('sor w = 1.3：')
[x, error_sor] = iterModel.sor_method(1.3);
[~, ~, loss_inf] = iterModel.cal_loss(x)
save_error_sor = [save_error_sor, error_sor];

% w = 1.8
rou_SOR = iterModel.cal_rou('S', 1.8)
disp('sor w = 1.8：')
[x, error_sor] = iterModel.sor_method(1.8);
[~, ~, loss_inf] = iterModel.cal_loss(x)
save_error_sor = [save_error_sor, error_sor];

%% draw different SOR
figure
lw = 2.5;
iter = 1:iterModel.iter_num;
plot(iter, save_error_sor(:,1), LineWidth=lw);
hold on
plot(iter, save_error_sor(:,2), LineWidth=lw);
plot(iter, save_error_sor(:,3), LineWidth=lw);
plot(iter, save_error_sor(:,4), LineWidth=lw);
ax = gca;
ax.LineWidth = 1.5;
ax.FontName = '微软雅黑';
ax.FontSize = 14;
grid on
set(gca,'GridLineStyle','--')
xlabel('迭代次数');
ylabel('误差');
title(['不同w的迭代法误差曲线','dim = ', num2str(dim)]);
legend('w=0.5', 'w=0.8', 'w=1.5', 'w=1.8');
