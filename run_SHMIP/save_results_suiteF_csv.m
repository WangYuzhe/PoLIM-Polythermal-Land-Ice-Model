% Date: 2019-5-29

clear,clc

% TIME SETTING
SPD = 86400; % [s]
SPY = 31536000; % [s]
dt = 1*SPD; % 0.5/24*SPD;
endTime = 1*SPY;
[arrayTime, numTimeStep] = set_time_step(dt, endTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_F1.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteF1_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteF1_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteF1_N.csv', suiteF1_N)
csvwrite('./results_csv/suiteF1_Qw.csv', suiteF1_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_F1 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_F1 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_F2.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteF2_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteF2_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteF2_N.csv', suiteF2_N)
csvwrite('./results_csv/suiteF2_Qw.csv', suiteF2_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_F2 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_F2 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_F3.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteF3_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteF3_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteF3_N.csv', suiteF3_N)
csvwrite('./results_csv/suiteF3_Qw.csv', suiteF3_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_F3 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_F3 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_F4.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteF4_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteF4_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteF4_N.csv', suiteF4_N)
csvwrite('./results_csv/suiteF4_Qw.csv', suiteF4_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_F4 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_F4 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_F5.mat
% arrayTime, upperN, upperN_std, middleN, middleN_std, lowerN, lowerN_std
suiteF5_N = [arrayTime,...
    mean(At_N(:,11:16),2)/1e6, std(At_N(:,11:16),0,2)/1e6,...
    mean(At_N(:,46:51),2)/1e6, std(At_N(:,46:51),0,2)/1e6,...
    mean(At_N(:,86:91),2)/1e6, std(At_N(:,86:91),0,2)/1e6];

% arrayTime, upperQw, middleQw, lowerQw
suiteF5_Qw = [arrayTime, At_Qw];

csvwrite('./results_csv/suiteF5_N.csv', suiteF5_N)
csvwrite('./results_csv/suiteF5_Qw.csv', suiteF5_Qw)

% time lag, amplitude of effective pressure
[~,time_max_rechg] = max(totalRechg);
[~,time_min_N] = min(mean(At_N,2));
timeLag_F5 = arrayTime(time_min_N) - arrayTime(time_max_rechg);
ampN_F5 = (max(mean(At_N,2)) - min(mean(At_N,2)))/1e6;

%===============================================================
timeLag_F = [timeLag_F1, timeLag_F2, timeLag_F3, timeLag_F4, timeLag_F5]/(24*3600*30);
ampN_F = [ampN_F1, ampN_F2, ampN_F3, ampN_F4, ampN_F5];

csvwrite('./results_csv/timeLag_F.csv',timeLag_F)
csvwrite('./results_csv/ampN_F.csv',ampN_F)