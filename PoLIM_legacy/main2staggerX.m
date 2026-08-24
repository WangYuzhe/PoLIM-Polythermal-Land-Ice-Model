function var_s = main2staggerX(var)
% Date: 2018-2-2

var_s = [var(:, 1), (var(:, 1:end-1) + var(:, 2:end))/2, var(:, end)];

end