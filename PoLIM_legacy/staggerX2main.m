function [var] = staggerX2main(var_s)
% Date: 2018-2-2

var = [var_s(:, 1), (var_s(:, 2:end-2)+var_s(:, 3:end-1))/2, var_s(:, end)];

end

