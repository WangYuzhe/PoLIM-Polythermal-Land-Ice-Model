function p = params()

%% physical constants
p.g = 9.81; % accel of gravity [m s-2]
p.R = 8.314; % universal gas constant [J mol K-1];
p.SPD = 86400; % day is this many seconds [s d-1]
p.SPY = 31556926; % year is this many seconds (i.e. 365.2422 days) [s yr-1]
p.rhoi = 910.0; % density of ice [kg m-3]
p.rhow = 1000.0; % water density [kg m-3]
p.n = 3.0; % exponent in Glen's flow law []
p.de0 = 1e-30; % small number in case of singularity [yr-2]
p.Lw = 3.34e5; % latent heat of fusion [J kg-1]
p.betaCC = 7.9e-8; % Clausius-Clapeyron constant [K Pa-1]
p.eta_w = 1.8e-3; % water viscosity (Hewitt&Schoof, 2017, 1.8e-3) [Pa s]
p.iter_max = 50; % maximum iterations for velocity solver []
p.Hmin = 0.1; % minimal ice thickness in case of singularity [m]
p.layers = 31; % vertical layers []

%% parameters related to thermodynamics
p.Tref = 223.15; % Reference temperature [K]
p.kc = 2.1; % cold ice conductivity [W m-1 K-1]
p.Cp = 2009; % ice specific heat capacity [J kg-1 K-1]
p.Kc = p.kc/(p.rhoi*p.Cp); % thermal diffusivity of cold ice [m2 s-1]
p.Kt = p.Kc*1e-3; % regularization of thermal diffusivity of temperate ice [m2 s-1]
p.k0 = 1e-12; % permeability factor (unconstrained, Hewitt&Schoof, 2017, Tab. 1) [m2]; range: 1e-12~5e-8
p.alpha_TEMP = 2.0; % exponent in compaction pressure model (Hewitt&Schoof, 2017); unconstrained, range: 2~3
p.Qgeo = 0.06; % geothermal heat flux [W m-2]
p.omega_max = 0.03; % upper bound of water content in temperate ice
p.is_greve_drain = 1; % water drainage in temperate ice using the Greve method
p.is_blatter_meltCTS = 1; % corrector step for cold layer (Blatter&Greve, 2015)

%% parameters related to Coulomb sliding law
p.lambda_max = 6; % wavelength of the dominant bedrock bumps [m]
p.m_max = 0.25; % maximum slope of the dominant bedrock bumps []; ref: 0.25
p.kflot = 0.1; % a fraction of flotation; N=Pi-Pw, Pw=kflot*Pi, N=(1-kflot)*Pi

%% parameters related to subglacial hydrology
p.is_subhydro = 0; % SUBGLACIAL HYDROLOGY
    % 1: couple to subglacial hydrology
    % 0: decouple to subglacial hydrology
p.lr = 2; % wavelength of bedrock bumps [m]; ref: 2
p.hr = 0.1; % height of bedrock bumps [m]; ref: 0.1
p.ktransm = 4e-4; % transmissivity coefficient [], larger value, larger N; ref: 4e-4
p.ev = 1e-5; % englacial void fraction (de Fleurian, 2018, eq. 10)
p.Hw_crit = 0.15; % critcal water sheet thickness (Pimentel & Flowers, 2010, Tab. 1)

%% options related to velocity solver
p.is_flowband = 1; % FLOWBAND MODE
    % true  (1): flowband mode
    % false (0): flowline mode

p.type_valley = 'trapz'; % VALLEY SHAPE
    % sves
    % ellipse
    % trapz (trapezoid)
    % rect (rectangle)

p.is_SBC_width = 0;
    % true  (1): incorporate flowband parameterization in stress-free condition
    % false (0): do not considier flowband parameterization in stress-free condition
    
p.f_nye = 1; % NYE SHAPE FACTOR []   

p.type_SBC = 'neum'; % SURFACE BOUNDARY CONDITION FOR FORCE BALANCE
    % neum: stress-free boundary condition (Neumann-type)
    % diri: specified surface velocity (Dirichlet-type)

p.type_BBC = 'zero'; % BASAL BOUNDARY CONDITION FOR FORCE BALANCE
    % zero: no-slip bed
    % CFlaw_polyT: Coulomb friction law, sliding depends on basal thermal state
    % CFlaw_isoT: Coulomb friction law, sliding everywhere
    % LFlaw: linear friction law
    % LFlaw_E2: linear friction law for ISMIP-HOM E2
    % LFlaw_simple: simplified linear friction law

p.type_LBC = 'zero'; % LATERAL BOUNDARY CONDITION
    % zero: zero velocity at the terminus
    % calv: calving front

%% options related to evolution solver
p.is_thk_evolv = 0; % ICE THICKNESS EVOLUTION MODE
    % true (1): ice thickness evolves
    % false (0): ice thickness is fixed

p.is_surf_relax = 0; % SURFACE RELAXATION MODE
    % true (1): surface relaxation
    % false (0): no surface relaxation

%% options related to thermodynamics solver
p.type_enth_solver = 'isoT'; % ENTHALPY SOLVER
    % SEGM: standard enthalpy gradient model (SEGM; Aschwanden et al., 2012)
    % MEGM: modified enthalpy gradient model (MEGM; Hewitt I. and Schoof C., 2017)
    % isoT: isothermal regime

p.is_auto_enth_BBC = 1;
    % true (1): basal thermal BC is automatically set depending on Aschwanden
    % decision chart (Aschwanden et al., 2012, JG, p451, fig. 5)
    % false (0): basal thermal BC is specified    
        
p.type_enth_BBC = 'COLD_base_dry'; % BASAL BOUNDARY CONDITION FOR ENERGY BALANCE
    % see Kleiner et al., 2015, TC, section 2.2
    % COLD_base_dry: Cold base (dry)
    % TEMP_base: Temperate base
    % TEMP_layer: Temperate ice layer at base
    % COLD_base_wet: Cold base (wet)

p.is_blatter_meltCTS = 0;
    % determine if the corrector step for the cold layer is applied
    % see Blatter-Greve (2015)
    % 1: Apply the corrector step
    % 0: Do not apply the corrector step

p.is_enth_trans = 1;
    % true (1): transient enthalpy solver
    % false (0): steady-state enthalpy solver
    
p.is_greve_drain = 1;
    % This option is only used for SEGM model.
    % see Aschwanden et al., 2012, JG, p450, fig. 4
    % 1: Use the Greve's drainage function
    % 0: Do not use the Greve's drainage function
    
%% options related to creep law
p.type_Arrhenius = 'cuffey';
    % greve: Greve (2009)
    % cuffey: Cuffey & Paterson (2010)
    % iverson: Schohn-Iverson et al. (2025)
    
p.is_duval = 1;
    % 1: considering the effect of water on temperate ice creep based on Duval-Lliboutry relation
    % 0: neglecting the effect of water on temperate ice creep

end