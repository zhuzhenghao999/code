function [predictedValues, MAPE] = predict_SFIGNGBM(data, Best_score, step)
    r = Best_score(1);
    u = Best_score(2);
    N = Best_score(3);


    A1 = data;
    
    n = length(A1);

    % 输入分数阶矩阵
    Ar = apply_fractional(r, A1);
    A2 = Ar;

    % 利用伯努利方法得出Y序列
    Y = A2.^(1 - N);

    % 构建数据矩阵
    for i = 2:n
        Yr(i) = Y(i) - Y(i - 1);
        Zy(i) = 0.5 * (Y(i) + Y(i - 1));
        R(i) = 0.5 * (gammainc(i, u) * gamma(u) + gammainc(i - 1, u) * gamma(u));
    end
    Zy(1) = [];
    Yr(1) = [];
    R(1) = [];
    B = (1 - N) .* [-Zy; R; ones(1, n - 1)];

    % 解出参数a,b,c
    C = pinv(B * B') * B * Yr';

    % 预测数据
    X_1 = predict_data(A2, C, Best_score, step);

    % 计算相对误差
    for i = 1:n
        APE(i) = abs(X_1(i) - A1(i)) / A1(i);
    end
    MAPE = sum(APE) / n;

    % 返回预测值和MAPE
    predictedValues = X_1;
end
