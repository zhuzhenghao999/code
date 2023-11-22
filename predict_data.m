function predicted_values = predict_data(af_data, C, best_parameters, step)
    a = C(1); b = C(2); c = C(3);
    A2 = af_data;
    n = length(A2);
    r = best_parameters(1);u = best_parameters(2) ;N =best_parameters(3);
    % 初始化变量

    gammaTermVector = zeros(1, n + step);
    expTermVector = zeros(1, n + step);
    totalMatrix = zeros(n + step, n + step);
    totalSumVector = zeros(1, n + step);
    totalExpVector = zeros(1, n + step);
    XrVector = zeros(1, n + step);
    
    % 拟合方程
    for j = 2:n + step
        for i = 2:j
            gammaTermVector(i) = gammainc(i, u) * gamma(u) + gammainc(i-1, u) * gamma(u) + 2 * c;
            expTermVector(i) = exp(-a * (1 - N) * (n + step - i + 0.5));
            totalMatrix(j, i) = gammaTermVector(i) * expTermVector(i);
        end
        totalSumVector(j) = sum(totalMatrix(j, :));
    end
    
    for i = 2:n + step
        totalExpVector(i) = exp(-a * (1 - N) * (i - 1)) * (A2(1)^(1 - N));
        XrVector(i) = (totalExpVector(i) + totalSumVector(i) * b * (1 - N) * 0.5)^(1 / (1 - N));
    end
        XrVector(1) = A2(1);

    % 分数阶还原数据
    predicted_values = inverse_fractional(r, XrVector);
end
