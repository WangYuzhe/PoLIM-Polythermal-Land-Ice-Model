function [arrayTime, numTimeStep] = set_time_step(dt, endTime)
% Date: 2019-10-21
% set time step
    % if dt = 2, endTime = 7
    % arrayTime: [1 3 5 7]
    % numTimeStep: 4
    % then arrayTime(3) = 5

arrayTime = dt:dt:(endTime);
arrayTime = arrayTime';
numTimeStep = length(arrayTime);

end