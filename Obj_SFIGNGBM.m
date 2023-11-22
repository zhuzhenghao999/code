function fitness = Obj_SFIGNGBM(x) 
% %预设数据阶数r，gamma函数参数u，幂指数nm
r =x(1);u =x(2) ;N =x(3);

%% 优化函数% 导入数据
%设置季节 k =1,2,3,4 (春夏秋冬）
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
a = C(1);b = C(2);c = C(3);
% 拟合方程
%g(i) = a0*(i-1)+a1*(i^(alpha+1)-1)/(alpha+1)
total = zeros(n,n);
for j = 2:n
    for i = 2:j
        to(i) = gammainc(i,u)*gamma(u)+gammainc(i-1,u)*gamma(u)+2*c;
        tal(i) =exp(-a*(1-N)*(n-i+0.5)) ;
        total(j,i) =(gammainc(i,u)*gamma(u)+gammainc(i-1,u)*gamma(u)+2*c)*exp(-a*(1-N)*(j-i+0.5));
    end
    total1(j) = sum(total(j,:));
end
for i = 2:n
    total2(i) = exp(-a*(1-N)*(i-1))*(A2(1)^(1-N));
    Xr(i) = (total2(i)+total1(i)*b*(1-N)*0.5).^(1/(1-N));
    Xr(1) = A1(1);
end
% 还原数据
X_1 = inverse_fractional(r,Xr);
for i = 1:n
    APE(i) = abs(X_1(i)-A1(i))/A1(i);
end
MAPE = sum(APE)/n;
fitness = MAPE;