function tests = test_smb_driver
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(project_root)
testCase.TestData.project_root = project_root;
end

function teardownOnce(testCase)
rmpath(testCase.TestData.project_root)
end

function testPddAnnualIntegrationAndRestart(testCase)
climate = synthetic_climate();
hS = [2600, 3000, 3400, 3800];
active = [true, true, true, false];
p = params();

driver = struct();
driver.mode = 'pdd';
driver.climate = climate;
driver.calendar_start_year = 2008;
driver.model_start_time = 1;
driver.hydro_start_month = 10;
driver.zref = 2943;
driver.rng_seed = 17;

rng_before = rng;
[smb1, state1] = run_smb(hS, active, driver, 1, p);
rng_after = rng;

testCase.verifyEqual(rng_after, rng_before)
testCase.verifySize(smb1, size(hS))
testCase.verifyTrue(all(isfinite(smb1)))
testCase.verifyEqual(smb1(~active), 0)
testCase.verifySize(state1.last.monthly.smb, [12, numel(hS)])
testCase.verifyEqual(state1.state.d_snow(~active), 0)

% A restart from the same driver state and seed must reproduce the year.
[smb_restart, state_restart] = run_smb(hS, active, driver, 1, p);
testCase.verifyEqual(smb_restart, smb1, 'AbsTol', 1e-12)
testCase.verifyEqual(state_restart.state.d_snow, state1.state.d_snow, ...
    'AbsTol', 1e-12)

% A repeated request is cached and must not advance snow or RNG state.
[smb_cached, state_cached] = run_smb(hS, active, state1, 1, p);
testCase.verifyEqual(smb_cached, smb1, 'AbsTol', 1e-12)
testCase.verifyEqual(state_cached.state, state1.state)

[smb2, state2] = run_smb(hS, active, state1, 2, p);
testCase.verifyTrue(all(isfinite(smb2)))
testCase.verifyEqual(state2.last.calendar_year, 2009)
end

function testMonthlySnowAndRainPartition(testCase)
zelev = 3000;
zref = 3000;
prcp = 1e-5;
n_days = 30;

[~, abl_cold, acc_cold, snow_cold] = calc_smb_pdd( ...
    -10, 0, prcp, 0, zelev, zref, 1, n_days);
testCase.verifyEqual(abl_cold, 0, 'AbsTol', 1e-12)
testCase.verifyGreaterThan(acc_cold, 0)
testCase.verifyEqual(snow_cold, acc_cold, 'AbsTol', 1e-12)

[~, ~, acc_warm, snow_warm] = calc_smb_pdd( ...
    10, 0, prcp, 0, zelev, zref, 7, 31);
testCase.verifyEqual(acc_warm, 0, 'AbsTol', 1e-12)
testCase.verifyEqual(snow_warm, 0, 'AbsTol', 1e-12)

[~, ~, acc_mixed] = calc_smb_pdd( ...
    1.5, 0, prcp, 0, zelev, zref, 4, n_days);
p_smb = params_smb();
prcp_month = prcp*n_days*86400/p_smb.rhow*p_smb.c_prcp;
testCase.verifyEqual(acc_mixed, prcp_month/2, 'RelTol', 1e-12)
end

function climate = synthetic_climate()
year = [2008*ones(3,1); 2009*ones(12,1); 2010*ones(9,1)];
month = [(10:12)'; (1:12)'; (1:9)'];
ndays = eomday(year, month);

% Smooth seasonal forcing in absolute temperature, matching the format in
% the original tests/test_run_smb.m example.
T_celsius = -5 - 10*cos(2*pi*(month-1)/12);
climate.date_glac = [year, month, ndays];
climate.tas_glac = T_celsius + 273.15;
climate.tas_std_glac = 2*ones(size(month));
climate.pr_glac = 1e-5*ones(size(month));
end
