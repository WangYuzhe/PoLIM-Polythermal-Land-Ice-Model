% Date: 2017-11-13
clc, clear
clearvars
clearvars -global

L_range = [5*1000, 10*1000, 20*1000, 40*1000, 80*1000, 160*1000];
n = length(L_range);

for i = 1:n
    L = L_range(i);
    func_main(L);
    
    clearvars -except L_range n i
    clearvars -global
end