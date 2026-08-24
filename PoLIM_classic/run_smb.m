function [smb, driver] = run_smb(hS, active, driver, model_time, p)
%RUN_SMB Calculate one year of surface mass balance for a classic flowline.
%
% The PDD mode integrates monthly climate in hydrological-year order and
% retains snow and random-number state in driver.state. In auto mode a
% missing climate file falls back to the configured legacy SMB method.
%
% Required inputs:
%   hS         surface elevation [m]
%   active     active main-grid mask
%   driver     configuration and persistent SMB state
%   model_time current model time [a]
%   p          PoLIM parameter structure
%
% Important driver fields:
%   mode                'auto', 'pdd', 'gradient', or 'fujita'
%   climate_file        MAT file containing gcm_glac (PDD/auto)
%   calendar_start_year climate year mapped to model_start_time
%   model_start_time    first model time [a]
%   hydro_start_month   hydrological-year start month (default 10)
%   zref                reference-station elevation [m]

hS = hS(:)';
active = logical(active(:)');
if numel(active) ~= numel(hS) || isempty(hS) || ...
        any(~isfinite(hS)) || ~isscalar(model_time) || ~isfinite(model_time)
    error('run_smb:InvalidState', ...
        'hS and active must be equal-length vectors and model_time finite.')
end
if ~isstruct(driver) || ~isstruct(p) || ~isfield(p, 'rhoi') || ...
        ~isscalar(p.rhoi) || ~isfinite(p.rhoi) || p.rhoi <= 0
    error('run_smb:InvalidConfiguration', ...
        'driver and p must be structures, and p.rhoi must be available.')
end

driver = apply_defaults(driver);

% Repeated requests at the same output time must not advance snow or RNG.
if isfield(driver, 'last') && isstruct(driver.last) && ...
        isfield(driver.last, 'model_time') && ...
        isequal(driver.last.model_time, model_time)
    smb = driver.last.smb;
    smb(~active) = 0;
    return
end

mode = lower(char(driver.mode));
if strcmp(mode, 'auto')
    [driver, has_climate] = load_climate(driver);
    if has_climate
        mode = 'pdd';
    else
        mode = lower(char(driver.fallback));
        if ~driver.warned_missing_climate
            warning('run_smb:MissingClimateFallback', ...
                ['PDD climate input was not found. Using the %s SMB ' ...
                 'fallback.'], mode)
            driver.warned_missing_climate = true;
        end
    end
elseif strcmp(mode, 'pdd')
    [driver, has_climate] = load_climate(driver);
    if ~has_climate
        error('run_smb:MissingClimate', ...
            'PDD mode requires driver.climate or a valid climate_file.')
    end
end

switch mode
    case 'pdd'
        [smb, driver] = integrate_pdd_year(hS, active, driver, model_time, p);
    case 'gradient'
        z_active = hS(active);
        if isempty(z_active)
            smb = zeros(size(hS));
        else
            smb = calc_smb_gradient(hS, median(z_active), driver.grad_smb);
        end
    case 'fujita'
        smb = calc_smb_fujita(driver.fujita_zref, driver.Tma, hS);
    otherwise
        error('run_smb:InvalidMode', ...
            'SMB mode must be auto, pdd, gradient, or fujita.')
end

smb = smb(:)';
smb(~active) = 0;
if any(~isfinite(smb))
    error('run_smb:NonfiniteSMB', 'The annual SMB contains nonfinite values.')
end

driver.last.model_time = model_time;
driver.last.smb = smb;
driver.last.mode = mode;

end

function driver = apply_defaults(driver)
defaults = struct( ...
    'mode', 'auto', ...
    'fallback', 'gradient', ...
    'climate_file', 'gfdl_hist.mat', ...
    'calendar_start_year', 2008, ...
    'model_start_time', 1, ...
    'hydro_start_month', 10, ...
    'zref', 2943, ...
    'fujita_zref', 2500, ...
    'Tma', -5.5, ...
    'grad_smb', 5e-2, ...
    'd_snow_init', 2e-2, ...
    'rng_seed', 0, ...
    'warned_missing_climate', false);

names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(driver, name) || isempty(driver.(name))
        driver.(name) = defaults.(name);
    end
end
if ~isfield(driver, 'state') || ~isstruct(driver.state)
    driver.state = struct();
end
if ~isfield(driver, 'last') || ~isstruct(driver.last)
    driver.last = struct();
end

valid_mode = ischar(driver.mode) || ...
    (isstring(driver.mode) && isscalar(driver.mode));
if ~valid_mode || ~isscalar(driver.hydro_start_month) || ...
        driver.hydro_start_month < 1 || driver.hydro_start_month > 12 || ...
        driver.hydro_start_month ~= fix(driver.hydro_start_month)
    error('run_smb:InvalidConfiguration', ...
        'The SMB mode and hydrological-year start month are invalid.')
end
end

function [driver, has_climate] = load_climate(driver)
has_climate = isfield(driver, 'climate') && isstruct(driver.climate) && ...
    ~isempty(fieldnames(driver.climate));
if has_climate
    return
end

has_file = (ischar(driver.climate_file) || ...
    (isstring(driver.climate_file) && isscalar(driver.climate_file))) && ...
    isfile(driver.climate_file);
if ~has_file
    return
end

data = load(driver.climate_file);
if isfield(data, 'gcm_glac') && isstruct(data.gcm_glac)
    driver.climate = data.gcm_glac;
else
    names = fieldnames(data);
    idx_struct = find(structfun(@isstruct, data), 1, 'first');
    if isempty(idx_struct)
        error('run_smb:InvalidClimateFile', ...
            'The climate MAT file does not contain a climate structure.')
    end
    driver.climate = data.(names{idx_struct});
end
has_climate = true;
end

function [smb, driver] = integrate_pdd_year(hS, active, driver, model_time, p)
climate = driver.climate;
required = {'date_glac', 'tas_glac', 'tas_std_glac', 'pr_glac'};
if ~all(isfield(climate, required))
    error('run_smb:InvalidClimate', ...
        'Climate must contain date_glac, tas_glac, tas_std_glac, and pr_glac.')
end

date = climate.date_glac;
if size(date, 2) < 2
    error('run_smb:InvalidClimate', ...
        'date_glac must contain at least year and month columns.')
end
year_all = date(:,1);
mon_all = date(:,2);
if size(date, 2) >= 3
    ndays_all = date(:,3);
else
    ndays_all = eomday(year_all, mon_all);
end
nclim = size(date, 1);
if numel(climate.tas_glac) ~= nclim || ...
        numel(climate.pr_glac) ~= nclim || ...
        ~ismember(numel(climate.tas_std_glac), [12, nclim]) || ...
        any(~isfinite(year_all)) || any(~isfinite(mon_all)) || ...
        any(mon_all < 1 | mon_all > 12 | mon_all ~= fix(mon_all)) || ...
        any(~isfinite(ndays_all) | ndays_all < 1 | ndays_all ~= fix(ndays_all))
    error('run_smb:InvalidClimate', ...
        'Climate arrays and calendar columns have inconsistent dimensions or values.')
end

year0 = driver.calendar_start_year + ...
    round(model_time - driver.model_start_time);
mon0 = driver.hydro_start_month;
if mon0 == 1
    idx = year_all == year0;
else
    idx = (year_all == year0 & mon_all >= mon0) | ...
        (year_all == year0 + 1 & mon_all < mon0);
end
idx = find(idx);
if numel(idx) ~= 12
    error('run_smb:MissingClimateYear', ...
        ['Expected 12 monthly records for the hydrological year starting ' ...
         'in %d-%02d, but found %d.'], year0, mon0, numel(idx))
end

[~, order] = sortrows([year_all(idx), mon_all(idx)], [1 2]);
idx = idx(order);
nz = numel(hS);
needs_snow = ~isfield(driver.state, 'd_snow') || ...
    numel(driver.state.d_snow) ~= nz;
if needs_snow
    driver.state.d_snow = zeros(1, nz);
    z_active = hS(active);
    if ~isempty(z_active)
        zmed = median(z_active);
        driver.state.d_snow(active & hS >= zmed) = driver.d_snow_init;
    end
else
    driver.state.d_snow = driver.state.d_snow(:)';
end
driver.state.d_snow(~active) = 0;

caller_rng = rng;
cleanup_rng = onCleanup(@() rng(caller_rng));
if isfield(driver.state, 'rng') && isstruct(driver.state.rng)
    rng(driver.state.rng)
else
    rng(driver.rng_seed, 'twister')
end

monthly = struct();
monthly.date = date(idx,:);
monthly.smb = zeros(12, nz);
monthly.abl = zeros(12, nz);
monthly.acc = zeros(12, nz);

for i = 1:12
    k = idx(i);
    Tair = climate.tas_glac(k);
    if Tair > 100 % test_run_smb stores absolute temperature [K]
        Tair = Tair - 273.15;
    end
    if numel(climate.tas_std_glac) == 12
        Tstd = climate.tas_std_glac(mon_all(k));
    else
        Tstd = climate.tas_std_glac(k);
    end
    prcp = climate.pr_glac(k);

    [mb, abl, acc, driver.state.d_snow] = calc_smb_pdd( ...
        Tair, Tstd, prcp, driver.state.d_snow, hS, driver.zref, ...
        mon_all(k), ndays_all(k));
    monthly.smb(i,:) = mb;
    monthly.abl(i,:) = abl;
    monthly.acc(i,:) = acc;
end
driver.state.rng = rng;
driver.state.d_snow(~active) = 0;

% calc_smb_pdd returns metres water equivalent; continuity evolves metres
% of ice, so convert the annual total with the density ratio.
p_smb = params_smb();
smb = sum(monthly.smb, 1)*(p_smb.rhow/p.rhoi);
smb(~active) = 0;
monthly.smb(:,~active) = 0;
monthly.abl(:,~active) = 0;
monthly.acc(:,~active) = 0;

driver.last.monthly = monthly;
driver.last.calendar_year = year0;

end
