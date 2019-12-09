% Date: 2019-11-10
function mouInput = calc_expB_mouInput(runID)
global M xi dx

num_mou_range = [1, 10, 20, 50, 100];
ice_L = 100e3;
ice_W = 20e3;
totalWaterInput = 4.5e-8*ice_L*ice_W; % [m3 s-1]
single_mouInput = totalWaterInput/num_mou_range(runID)/(dx*ice_W);

% load moulin information
csvFileName = strcat('B', num2str(runID), '_M.csv');
dataMou = csvread(csvFileName);
coord_x = ice_L - dataMou(:,2);

mouInput = zeros(1,M);
for i = 1:M
    index = find(coord_x == xi(i));
    if isempty(index)
        mouInput(1,i) = 0;
    else
        mouInput(1,i) = length(index)*single_mouInput;
    end
end
end