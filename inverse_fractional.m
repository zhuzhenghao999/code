function predicted_value = inverse_fractional(r, value)
    if r == 1
        predicted_value = value;
    else
        n = length(value);
        processed_value = zeros(size(value));

        for i = 1:n
            for j = 1:n
                tmp = gamma(i-j+r-1+1)/(gamma(r-1+1)*gamma(i-j+1));
                processed_value(i) = processed_value(i) + tmp * value(j);
            end
        end

        predicted_value = processed_value - circshift(processed_value, [0, 1]);
        predicted_value(1) = value(1);
    end
end
