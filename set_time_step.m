function [arrayTime, numTimeStep] = set_time_step(dt,startTime,endTime)
% Date: 2019-10-21
% set time step
    % if dt = 2, startTime = 2010, endTime = 2020
    % then arrayTime = [2010,2012,2014,2016,2018,2020], numTimeStep: 6
    % then arrayTime(2) = 2012

arrayTime = startTime:dt:endTime;
arrayTime = arrayTime';
numTimeStep = length(arrayTime);

end