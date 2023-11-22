clear
clc
rng('shuffle')
%%  参数设置
pop = 300;              % 种群数目
Max_iter = 30;         % 迭代次数
dim = 3;               % 优化参数个数
lb = [0,   0,  0];       % 下限
ub = [10, 10,  10];       % 上限

fobj = @(x)Obj_SFIGNGBM(x) ;

%% 利用鲸鱼算法求出最优参数 r u N
[Best_pos, Best_score, curve] = WOA(pop, Max_iter, lb, ub, dim, fobj);
%% 得出模型最优参数进行预测
r = Best_score(1);u = Best_score(2) ;N =Best_score(3);

%% 导入数据

%设置季节 k =1,2,3,4 (春夏秋冬）（以及Obj_SFIGNGBM()季节）
k = 2;
A = readmatrix("数据0129.xlsx",'Sheet','USA','Range','B2:B31');
for i = 1:8
    j = 4*i-k;
    A1(i) = A(j);
end
% A1 = A;
n = length(A1);
% 输入分数阶矩阵
Ar = apply_fractional(r,A1);
A2 = Ar;

% 利用伯努利方法得出Y序列
Y = A2.^(1-N);
% 构建数据矩阵
for i = 2:n
    Yr(i) = Y(i)-Y(i-1);
    Zy(i) = 0.5*(Y(i)+Y(i-1));
    R(i) = 0.5*(gammainc(i,u)*gamma(u) +gammainc(i-1,u)*gamma(u));
end
Zy(1) = [];Yr(1) = [];R(1) = [];
B = (1-N).*[-Zy;R;ones(1,n-1)];
% 解出参数a,b,c
C = pinv(B*B')*B*Yr';
% 设置预测步长
step = 2;
% 预测
X_1 = predict_data(A2, C, Best_score, step);
for i = 1:n
APE(i) = abs(X_1(i)-A1(i))/A1(i);
end
MAPE = sum(APE)/n;
fitness = MAPE;

% 打印输出 MAPE
fprintf('MAPE: %.4f%%\n', MAPE*100);
season = ['春', '夏', '秋', '冬'];
% 打印输出季节信息
fprintf('季节: %s\n', season(k));
% 打印输出参数信息
fprintf('参数分数阶阶参数r: %.4f\n', r);
fprintf('参数gamma函数参数u: %.4f\n', u);
fprintf('伯努利参数N: %.4f\n', N);
%% 做图
load('color_list')
curve_compare=[];     %算法的过程函数比较
Fival_compare=[];  %算法的最终目标比较
curve_compare=[curve_compare;curve];
Fival_compare=[Fival_compare,Best_pos];
load('color_list')
figure(2)
color_all=color_list(randperm(length(color_list)),:);

%画迭代过程曲线
for N=1:length(Fival_compare)
     plot(curve_compare(N,:),'Color',color_all(N,:),'LineWidth',2)
     hold on
end
xlabel('迭代次数');
ylabel('目标函数值');
grid on
box on
