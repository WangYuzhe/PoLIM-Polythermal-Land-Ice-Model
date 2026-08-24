function p = params_smb()

p.f_snow = 2e-3; % 1.75~4.5 [m d-1 K-1], Zekollari et al (2019), p1129
p.c_prcp = 1.8; % Huss&Hock (2015), p8; range: 0.8~2.0, Zekollari et al (2019), p1129
p.grad_prcp = 1.5e-4; % [m] 1~2.5%/100 m, 1e-4~2.5e-4 (Huss2015, eq.2)
p.Tmelt = 0; % [Celsius] temperature when glacier melt occurs
p.Tsl = 1.5; % [Celsius] temperature threshold for precipitation differentiation
p.rhow = 1000; % [kg m-3]

p.SMB_grad = 4e-3;
p.SMB_max = 1.5;

end