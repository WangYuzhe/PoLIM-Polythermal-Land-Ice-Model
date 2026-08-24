clear, clc

% predefined parameters
indir = 'E:\Mix\WYZ_workspace\thermal_modeling\manu_lhg12_enthalpy\data';
outdir = './geo_inputs';
glac_name = 'lhg12_glate';

% initialized parameters
is_smooth = 1;
is_width = 1;
is_uobs = 0;

% derived paths
file_attrs_geom = dir(fullfile(indir, 'geom*.csv'));
incsv_geom = fullfile(indir, file_attrs_geom.name);
outmat_geo = fullfile(outdir, strcat('geo_', glac_name,'.mat'));

if is_uobs
    file_attrs_uobs = dir(fullfile(indir, 'uobs*.csv'));
    incsv_uobs = fullfile(indir, file_attrs_uobs.name);
    uobs = csvread(incsv_uobs,1,1);
end

% read csv files
geom = csvread(incsv_geom,1,1);
geo.xi = geom(1,:);

% take the original surface and bedrock elevations
geo.hB = geom(2,:);
geo.hS = geom(3,:);
geo.H = geom(4,:);

if is_width
    geo.Wsurf = geom(8,:);
end

% take the smooth surface and bedrock elevations
if is_smooth
    geo.hB = geom(5,:);
    geo.hS = geom(6,:);
    geo.H = geom(7,:);
end

% save uobs
if is_uobs
    save(outmat_geo, 'geo', 'uobs');
else
    save(outmat_geo, 'geo');
end


figure(1)
plot(geo.xi, geo.hB, 'k-')
hold on
plot(geo.xi, geo.hS, 'b-')