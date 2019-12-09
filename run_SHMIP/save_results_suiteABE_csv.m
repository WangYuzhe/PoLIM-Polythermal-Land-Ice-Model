% Date: 2019-5-19

clear,clc

% TOPOGRAPHY
xi_icesheet = 0:1000:100e3;
xi_glacier = 0:60:6e3;

%=========================Suite A=========================
load N_A1.mat
load N_A2.mat
load N_A3.mat
load N_A4.mat
load N_A5.mat
load N_A6.mat

suiteA_N = [xi_icesheet; N_A1; N_A2; N_A3; N_A4; N_A5; N_A6];
suiteA_N = suiteA_N';
csvwrite('./results_csv/suiteA_N.csv', suiteA_N)

clearvars -EXCEPT xi_icesheet xi_glacier

%=========================Suite B=========================
load N_B1.mat
load N_B2.mat
load N_B3.mat
load N_B4.mat
load N_B5.mat

suiteB_N = [xi_icesheet; N_B1; N_B2; N_B3; N_B4; N_B5];
suiteB_N = suiteB_N';
csvwrite('./results_csv/suiteB_N.csv', suiteB_N)

clearvars -EXCEPT xi_icesheet xi_glacier

%=========================Suite E=========================
load N_E1.mat
load N_E2.mat
load N_E3.mat
load N_E4.mat
load N_E5.mat

suiteE_N = [xi_glacier; N_E1; N_E2; N_E3; N_E4; N_E5];
suiteE_N = suiteE_N';
csvwrite('./results_csv/suiteE_N.csv', suiteE_N)

clearvars -EXCEPT xi_icesheet xi_glacier