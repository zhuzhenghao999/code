clc
clear
% rng('shuffle')
%%  参数设置
pop = 800;              % 种群数目
Max_iter = 30;         % 迭代次数
dim = 3;               % 优化参数个数
lb = [0,   0,  0];       % 下限
ub = [10, 10,  10];       % 上限
load('data.mat');
step = 2;
groupSize = 4;
data;
padding = mod(length(data), groupSize);
data = vertcat(data, NaN(padding, 1));
groupedData = reshape(data, groupSize, [])';

groups = cell(1, size(groupedData, 2));
results = struct('predictedValues', {}, 'MAPE', [], 'Best_score', []);
% 循环处理每个分组
for i = 1:size(groupedData, 2)
    % 从 groupedData 中提取当前分组
    currentGroup = groupedData(:, i);
    

    groups{i} = currentGroup(~isnan(currentGroup));
    

%%
groupdata = groups{i};
fobj = @(x)Obj_SFIGNGBM(x,groupdata) ;



[Best_pos, Best_score, curve] = WOA(pop, Max_iter, lb, ub, dim, fobj);


%% 得出模型最优参数进行预测

[predictedValues, MAPE] = predict_SFIGNGBM(groupdata, Best_score, step);

    numRetry = 0;
    while any(predictedValues < 0) && numRetry < 5
        disp('Detected negative values in predictedValues. Re-running optimization...'); 
        [Best_pos, Best_score, curve] = WOA(pop, Max_iter, lb, ub, dim, fobj);       
        [predictedValues, MAPE] = predict_SFIGNGBM(groupdata, Best_score, step);
        numRetry = numRetry + 1;
    end
    % 将结果保存到结构体数组中
    results(i).predictedValues = predictedValues;
    results(i).MAPE = MAPE;
    results(i).Best_score = Best_score;
    %% 做图
    load('color_list')
    curve_compare=[];     %算法的过程函数比较
    Fival_compare=[];  %算法的最终目标比较
    curve_compare=[curve_compare;curve];
    Fival_compare=[Fival_compare,Best_pos];
    load('color_list')
    figure(2)
    color_all=color_list(randperm(length(color_list)),:);
    
    % 画迭代过程曲线
    for N=1:length(Fival_compare)
         plot(curve_compare(N,:),'Color',color_all(N,:),'LineWidth',2)
         hold on
    end
    xlabel('迭代次数');
    ylabel('目标函数值');
    grid on
    box on
end

% 结果的访问示例
for i = 1:size(results, 2)
    disp(['Group ', num2str(i)]);
    disp('Predicted Values:'); 
    disp(results(i).predictedValues);
    disp(['MAPE: ', num2str(results(i).MAPE)]);
    disp('Best Score:');
    disp(results(i).Best_score);
end

% 初始化 data_2
data_2 = zeros(1, 4 * length(predictedValues) + step);

% 将预测值整合到 data_2 中
for i = 1:groupSize
    index = i:4:(i + 4 * (length(results(i).predictedValues) - 1) + step);
    data_2(index) = results(i).predictedValues;
end
disp(data_2);



% 原始数据
originalData = data;

% 收集所有分组的 MAPE
mapeValues = zeros(1, groupSize);
for i = 1:groupSize
    mapeValues(i) = results(i).MAPE;
end

% 计算平均 MAPE
averageMape = mean(mapeValues)*100;

% 绘制原始数据和预测数据
figure;
plot(originalData, '-o', 'LineWidth', 2, 'DisplayName', 'Original Data');
hold on;
plot(data_2, '-s', 'LineWidth', 2, 'DisplayName', 'Predicted Data');
hold off;

% 添加图例和标题
legend('Location', 'best');
title(['Comparison of Original and Predicted Data (Average MAPE: ', num2str(averageMape), '%)']);
xlabel('Time');
ylabel('Value');
grid on;
