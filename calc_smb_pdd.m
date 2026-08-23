function [m_b,m_abl,m_acc,out_d_snow] = calc_smb_pdd(Tair_ref,Tair_std_ref,prcp_ref,in_d_snow,zelev,zref,ith_mon,n_days)

if ~isscalar(Tair_ref) || ~isfinite(Tair_ref) || ...
        ~isscalar(Tair_std_ref) || ~isfinite(Tair_std_ref) || ...
        Tair_std_ref < 0 || ~isscalar(prcp_ref) || ...
        ~isfinite(prcp_ref) || prcp_ref < 0
    error('calc_smb_pdd:InvalidClimate', ...
        'Temperature, variability, and precipitation inputs must be finite scalars.')
end
if ~isscalar(ith_mon) || ith_mon < 1 || ith_mon > 12 || ...
        ith_mon ~= fix(ith_mon) || ~isscalar(n_days) || n_days < 1 || ...
        n_days ~= fix(n_days)
    error('calc_smb_pdd:InvalidCalendar', ...
        'ith_mon and n_days must be valid positive calendar integers.')
end

zelev = zelev(:)';
in_d_snow = in_d_snow(:)';
if numel(in_d_snow) ~= numel(zelev) || ...
        any(~isfinite(zelev)) || any(~isfinite(in_d_snow)) || ...
        any(in_d_snow < 0) || ~isscalar(zref) || ~isfinite(zref)
    error('calc_smb_pdd:InvalidGeometryState', ...
        'Snow state and glacier elevations must be finite vectors of equal size.')
end

%% physical parameters
%
p = params_smb();
f_snow = p.f_snow; % degree-day factor of snow
c_prcp = p.c_prcp; % scaling factor (Huss2015, eq.2)
grad_prcp = p.grad_prcp; %
Tmelt = p.Tmelt; % [Celsius]
Tsl = p.Tsl; % [Celsius]
rhow = p.rhow; % [kg m-3]
lapse_rate_mon = -6.5e-3*ones(1,12);

% derived parameters
f_ice = 2*f_snow;
n_xgrid = length(zelev);
Tsl_snow = Tsl - 1;
Tsl_rain = Tsl + 1;
total_secs = n_days*86400; % kg m-2 month-1

%%
% get the distributed monthly air temperature
Tair_glac = Tair_ref + (zelev-zref)*lapse_rate_mon(ith_mon); % <1*n_xgrid>

% Draw daily temperature anomalies about the monthly mean. The SMB driver
% owns and restores the random state, making this stochastic forcing
% reproducible across model restarts.
Tair_daily_ref = Tair_ref + Tair_std_ref.*randn(n_days,1);

Tair_daily_glac = Tair_daily_ref*ones(1,n_xgrid) +...
    ones(n_days,1)*(zelev - zref).*lapse_rate_mon(ith_mon); %<n_days*n_xgrid>

% get the distributed monthly precipitation
prcp_ref = prcp_ref*total_secs/rhow; % [m month-1]
prcp_glac = prcp_ref*c_prcp*max(1 + (zelev-zref)*grad_prcp, 0); % <1*n_xgrid>

% differentiate between solid and liquid precipitation
prcp_glac_snow = zeros(1,n_xgrid);
for i=1:n_xgrid
    if Tair_glac(i)<=Tsl_snow
        prcp_glac_snow(i) = prcp_glac(i);
    elseif Tair_glac(i)>=Tsl_rain
        prcp_glac_snow(i) = 0;
    else
        rain_fraction = (Tair_glac(i)-Tsl_snow)/...
            (Tsl_rain-Tsl_snow);
        prcp_glac_snow(i) = (1-rain_fraction)*prcp_glac(i);
    end
end

prcp_daily_glac_snow = prcp_glac_snow/n_days; % <1*n_xgrid>
%% CALCULATION
% calculate the monthly ablation using degree-day method along the flowline
logic_melt = (Tair_daily_glac>Tmelt); % boolean matrix, <n_days*n_xgrid>
pdd = logic_melt.*(Tair_daily_glac - Tmelt); % positive degree-days, <n_days*n_xgrid>

% potential PDD for melting away the whole snow
m_abl_daily = zeros(n_days,n_xgrid);
for i=1:n_days
    pot_pdd = in_d_snow/f_snow;
    for j=1:n_xgrid
        if in_d_snow(j)>0
            if pdd(i,j)<pot_pdd(j)
                m_abl_daily(i,j) = f_snow*pdd(i,j);
                in_d_snow(j) = in_d_snow(j) - m_abl_daily(i,j);
            else
                m_abl_daily(i,j) = f_snow*pot_pdd(j) +...
                    f_ice*(pdd(i,j) - pot_pdd(j));
                in_d_snow(j) = 0;
            end
        else
            m_abl_daily(i,j) = f_ice*pdd(i,j);
        end
    end
    in_d_snow(in_d_snow<0) = 0; % ensure that snow depth is not negative
    in_d_snow = in_d_snow + prcp_daily_glac_snow;
end

% calculate the monthly mass balance
m_abl = sum(m_abl_daily);
m_acc = prcp_glac_snow;
m_b = m_acc - m_abl;

% Snowfall has already been accumulated at daily resolution.
out_d_snow = in_d_snow;
end
