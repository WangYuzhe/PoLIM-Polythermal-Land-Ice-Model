function drainageRate = drainageFunc(omega)
% Date: 2019/4/19
% Author: Wang Yuzhe

if omega <= 0.01
    drainageRate = 0;
elseif omega > 0.01 && omega <=0.02
    drainageRate = 0.5*omega - 0.005;
elseif omega > 0.02 && omega <= 0.03
    drainageRate = 4.5*omega - 0.085;
else
    drainageRate = 0.05;
end

end
