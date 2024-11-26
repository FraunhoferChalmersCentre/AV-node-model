function result = logsum(values)

% Calculates the sum of exponentials logarithmically 

max_value = max(values);
sum_exp = sum(exp(values - max_value));
result = max_value + log(sum_exp);

end
