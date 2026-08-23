clear, clc

load('./geo_inputs/geo_yulong.mat');
load('gfdl_hist.mat');

xi = geo.xi;
zelev = geo.hS;
zref = 2943;
zmed = median(zelev);
lapse_rate_mon = -6.5e-3*ones(1,12);
in_d_snow = zeros(1,length(zelev));
in_d_snow(zelev>zmed) = 1e-2; % 5 cm

start_year = 2008;
end_year = 2014;

n_xgrid = length(zelev);
%% meteorological inputs
%
date_glac = gcm_glac.date_glac;
gcm_year = date_glac(:,1);
gcm_mon = date_glac(:,2);

% hydrological years
idx_start_year = find((gcm_year==start_year & gcm_mon==10));
idx_end_year = find((gcm_year==end_year & gcm_mon==9));

% subset the GCM input
Tair_ref = gcm_glac.tas_glac(idx_start_year:idx_end_year) - 273.15; % [Celsius], <n_mon>
Tair_std_ref = gcm_glac.tas_std_glac(idx_start_year:idx_end_year); % [K], <n_mon>
prcp_ref = gcm_glac.pr_glac(idx_start_year:idx_end_year); % [kg m-2 s-1], <n_mon>

gcm_mon = date_glac(idx_start_year:idx_end_year,2); % <n_mon>
gcm_n_days_in_month = date_glac(idx_start_year:idx_end_year,3); % <n_mon>

n_mon = length(gcm_mon);
%%
m_b_period = zeros(n_mon,n_xgrid);
m_abl_period = zeros(n_mon,n_xgrid);
m_acc_period = zeros(n_mon,n_xgrid);

for i=1:n_mon
    n_days = gcm_n_days_in_month(i);
    ith_mon = gcm_mon(i);
    
    if i==1
        in_d_snow = zeros(1,length(zelev));
        in_d_snow(zelev>zmed) = 2e-2; % 5 cm
    else
        in_d_snow = out_d_snow;
    end    
    [m_b,m_abl,m_acc,out_d_snow] = calc_smb_pdd(Tair_ref(i),Tair_std_ref(i),prcp_ref(i),in_d_snow,zelev,zref,ith_mon,n_days);
    m_b_period(i,:) = m_b;
    m_abl_period(i,:) = m_abl;
    m_acc_period(i,:) = m_acc;
end

figure
subplot(1,2,1)
plot(zelev, sum(m_b_period))

subplot(1,2,2)
plot(xi, sum(m_b_period))
