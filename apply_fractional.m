%%分数阶
function p_value = apply_fractional(r,value)
p_value = zeros(length(value));
if r == 1
    p_value = (cumsum(value))';
elseif r == 0
    p_value = value';
else
    for i = 1:length(value)
        for j = 1:length(value)
            tmp = gamma(i-j+r-1+1)/(gamma(r-1+1)*gamma(i-j+1));
            p_value(i) = p_value(i)+tmp*value(j);
        end
    end
end
