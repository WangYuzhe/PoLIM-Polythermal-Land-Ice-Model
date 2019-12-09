% 2019-6-17

clear, clc

ice_L = 100e3;
dx = 1000;
xi = 0:dx:ice_L;

load tune_result_N_A3.mat
ref_tune = xlsread('ref_tune.xlsx', 'A3');
ref_tune = ref_tune';


[n_row, n_col] = size(tune_result_N);
norm_delta_N = zeros(n_row, n_col);
interp_Nij = zeros(n_row*n_col, size(ref_tune,2));
for i = 1:n_row
    for j = 1:n_col
        Nij = tune_result_N{i,j};
        
        Nij = fliplr(Nij);
        
        Nij_1 = interp1(xi, Nij, ref_tune(1,:));
        
        interp_Nij((i-1)*n_row + j, :) = Nij_1;
        
        delta_N = Nij_1 - ref_tune(2,:);
        
        %std_delta_N(i,j) = std(delta_N);
        
        norm_delta_N(i,j) = norm(delta_N);
    end
end
[~,index_row] = min(norm_delta_N);
[~,index_col] = min(min(norm_delta_N));

best_row = index_row(index_col);
best_col = index_col;

k0_range = [1e-6, 1e-5, 1e-4, 2e-4, 4e-4, 5e-4, 1e-3, 1e-2, 1e-1];
etai_range = [1e12, 2e12, 3e12, 4e12, 1e13, 5e13, 1e14, 5e14, 1e15];

fprintf('The best combination of k0 and etai is (%3.1e, %3.1e) \n', k0_range(best_row), etai_range(best_col))

index_best_interp_Nij = (best_row-1)*n_row + best_col;

figure
hold on
FL=zeros(2,1);
FL(1) = plot(ref_tune(1,:)/1e3, ref_tune(2,:)/1e6, 'k');
FL(2) = plot(ref_tune(1,:)/1e3, interp_Nij(index_best_interp_Nij,:)/1e6, 'r');
legend(FL, 'Ref', 'PoLIM')
box on

xlabel('x (km)')
ylabel('N (MPa)')