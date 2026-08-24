% Date: 2019-5-29

clear,clc

% TIME SETTING
SPD = 24*3600; % [s d-1]
SPY = 365*SPD; % [s yr-1]
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_D1.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteD1_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteD1_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteD1_N.csv', suiteD1_N)
csvwrite('./results_csv/suiteD1_Qw.csv', suiteD1_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_D1 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_D1 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_D2.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteD2_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteD2_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteD2_N.csv', suiteD2_N)
csvwrite('./results_csv/suiteD2_Qw.csv', suiteD2_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_D2 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_D2 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_D3.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteD3_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteD3_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteD3_N.csv', suiteD3_N)
csvwrite('./results_csv/suiteD3_Qw.csv', suiteD3_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_D3 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_D3 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_D4.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteD4_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteD4_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteD4_N.csv', suiteD4_N)
csvwrite('./results_csv/suiteD4_Qw.csv', suiteD4_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_D4 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_D4 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_D5.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteD5_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteD5_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteD5_N.csv', suiteD5_N)
csvwrite('./results_csv/suiteD5_Qw.csv', suiteD5_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_D5 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_D5 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%===============================================================
timeLag_D = [timeLag_D1, timeLag_D2, timeLag_D3, timeLag_D4, timeLag_D5]/(24*3600*30);
ampN_D = [ampN_D1, ampN_D2, ampN_D3, ampN_D4, ampN_D5];

csvwrite('./results_csv/timeLag_D.csv', timeLag_D)
csvwrite('./results_csv/ampN_D.csv',ampN_D)