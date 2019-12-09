function [LT, RT] = basalBC_temperate_layer(LT1, LT2, LT3, RT)
global N

LT1(1) = 0;
LT2(1) = -2;
LT3(1) = 3;
RT(1) = 0;
LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);
LT(1,3) = -1;

% LT1(1) = 0;
% LT2(1) = -1/dzeta;
% LT3(1) = 1/dzeta;
% RT(1) = 0;
% LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);


end