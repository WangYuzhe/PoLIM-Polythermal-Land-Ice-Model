
global xi hB hS H

geo0 = load('./geo_inputs/geo_kangxiong1.mat');
geo0 = geo0.geo;

geo.xi = xi;
geo.hB = hB;
geo.hS = hS;
geo.H = H;
geo.Wsurf = geo0.Wsurf';

save('./geo_inputs/geo_kangxiong_relax.mat', 'geo')

