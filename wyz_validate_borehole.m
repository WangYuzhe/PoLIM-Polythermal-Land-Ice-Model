
global zeta

file_borehole_T = './geo_inputs/LaohugouNo12_borehole_Temp.csv';
nodes = [10, 11];
xi_core = 1819.5558067798636;

T_node1 = T(:, 10) - 273.15;
T_node2 = T(:, 11) - 273.15;
H_node1 = H(10);
H_node2 = H(11);
z_node1 = H_node1*(1 - zeta);
z_node2 = H_node2*(1 - zeta);


d1 = xi_core - xi(10);
d2 = xi(11) - xi_core;
T_bh = (T_node1*d2 + T_node2*d1) / (d1 + d2);

H_bh = (H_node1*d2 + H_node2*d1) / (d1 + d2);
z_bh = H_bh*(1 - zeta);

data_borehole = readmatrix(file_borehole_T);

figure
hold on
plot(data_borehole(:,2), data_borehole(:,1), 'k-o')
plot(T_bh, z_bh, 'r-o')
plot(T_node1, z_node1, 'b-')
plot(T_node2, z_node2, 'g-')
set(gca, 'YDir', 'reverse')




