classdef plot_benchmark
    methods(Static)
        function plot_shmip_D(timeSpan, At_N)
            xx = [timeSpan; flipud(timeSpan)];
            yy_high = [mean(At_N(:,11:16),2)/1e6+std(At_N(:,11:16),0,2)/1e6;...
                flipud(mean(At_N(:,11:16),2)/1e6-std(At_N(:,11:16),0,2)/1e6)];
            yy_middle = [mean(At_N(:,46:51),2)/1e6+std(At_N(:,46:51),0,2)/1e6;...
                flipud(mean(At_N(:,46:51),2)/1e6-std(At_N(:,46:51),0,2)/1e6)];
            yy_low = [mean(At_N(:,86:91),2)/1e6+std(At_N(:,86:91),0,2)/1e6;...
                flipud(mean(At_N(:,86:91),2)/1e6-std(At_N(:,86:91),0,2)/1e6)];
            
            hold on
            fill(xx, yy_high, 'r', 'Linestyle', 'none', 'FaceAlpha',0.5)
            fill(xx, yy_middle, 'c', 'Linestyle', 'none', 'FaceAlpha',0.5)
            fill(xx, yy_low, 'b', 'Linestyle', 'none', 'FaceAlpha',0.5)
            plot(timeSpan, mean(At_N(:,11:16),2)/1e6, 'r-', 'LineWidth', 1)
            plot(timeSpan, mean(At_N(:,46:51),2)/1e6, 'c-', 'LineWidth', 1)
            plot(timeSpan, mean(At_N(:,86:91),2)/1e6, 'b-', 'LineWidth', 1)
            
            box on
            
            xlim([90, 365])
            ylim([-5 10])
            xlabel('Month')
            ylabel('N (MPa)')
            set(gca, 'Xtick', [120, 181, 243, 304, 365])
            set(gca, 'XtickLabel', [4, 6, 8, 10, 12])
            set(gca, 'Ytick', [-4,0,4,8])
            set(gca, 'YtickLabel', [-4,0,4,8])
            
        end
        
        function plot_shadedErrorBar(data_NFS, data_FS, ylabel_text)
            % Date: 2017-9-13
            % Author: Wang Yuzhe
            % FS: full-Stokes
            % NFS: non-full-Stokes
            
            % NFS data
            x_NFS = data_NFS(:,1); % horizontal distance
            var_NFS = data_NFS(:,2); % variable (e.g., surface velocity, basal stress)
            var_std_NFS = data_NFS(:,3); % standard deviation
            
            var_upper_NFS = var_NFS + var_std_NFS;
            var_lower_NFS = var_NFS - var_std_NFS;
            
            % FS data
            x_FS = data_FS(:,1);
            var_FS = data_FS(:,2);
            var_std_FS = data_FS(:,3); % standard deviation
            
            var_upper_FS = var_FS + var_std_FS;
            var_lower_FS = var_FS - var_std_FS;
            
            
            F = zeros(4,1);
            hold on
            F(2) = fill([x_NFS; flipud(x_NFS)], [var_upper_NFS; flipud(var_lower_NFS)],...
                [0.929, 0.694, 0.125], 'linestyle', 'none', 'FaceAlpha', 0.5);
            F(4) = plot(x_NFS, var_NFS, 'color', [0.85, 0.325, 0.098], 'linewidth', 1);
            
            F(1) = fill([x_FS; flipud(x_FS)], [var_upper_FS; flipud(var_lower_FS)],...
                [0.301, 0.745, 0.933], 'linestyle', 'none', 'FaceAlpha', 0.5);
            F(3) = plot(x_FS, var_FS, 'color', [0, 0.447, 0.741], 'linewidth', 1);
            
            % F(1).FaceAlpha = 0.5;
            % alpha(F(1), 0.7)
            % alpha(F(3), 0.7)
            
            hold off
            
            legend(F, 'FS', 'NFS', 'FS Mean', 'NFS Mean', 'location', 'northwest')
            
            xlabel('Normalized x')
            ylabel(ylabel_text)
            
            box on
            grid on
        end
        
        function plot_enthalpy_ExpA(arrayTime,At_T,At_m_basal,At_thk_w,At_isTEMP,enthA_analy_result)
            subplot(2,2,1)
            plot(arrayTime/1000, At_T(1,:)-273.15, 'k', 'linewidth', 1)
            xlim([0 310])
            ylim([-30 5])
            ylabel('T_b (^\circC)')
            grid on
            box on
            
            subplot(2,2,2)
            hold on
            plot(enthA_analy_result(:,1)/1000, enthA_analy_result(:,2), 'r', 'linewidth', 2)
            plot(arrayTime/1000, At_m_basal(:,1)*1000, 'k', 'linewidth', 1)
            xlim([0 310])
            ylim([-4 4])
            ylabel('a_b (mm a^{-1} w.e.)')
            grid on
            box on
            legend('Analytical', 'PoLIM', 'location', 'southwest')
            
            subplot(2,2,3)
            plot(arrayTime/1000, At_thk_w(:,1), 'k', 'linewidth', 1)
            xlim([0 310])
            ylim([-10 160])
            ylabel('H_w (m)')
            xlabel('Time (ka)')
            grid on
            box on
            
            subplot(2,2,4)
            plot(arrayTime/1000, At_isTEMP(:,1), 'k', 'linewidth', 1)
            xlim([0 310])
            ylim([-0.2 1.2])
            ylabel('is temperate base')
            xlabel('Time (ka)')
            grid on
            box on
            
            set(gcf, 'position', [246 114 831 536])
        end
        
        function plot_enthalpy_ExpB(zeta,E,T,omega,enthB_analy_z,...
                enthB_analy_E,enthB_analy_T,enthB_analy_omega)
            subplot(1,3,1)
            hold on
            plot(enthB_analy_E/1000, enthB_analy_z, 'r', 'linewidth', 3)
            plot(E(:,1)/1000, zeta, '-k', 'linewidth', 1.5)
            % plot(E(jcts,:,end)/10^3, zeta(jcts), 'ob')
            xlim([94,108])
            set(gca, 'xtick', [96, 100, 104, 108])
            set(gca, 'ytick', 0:0.1:1)
            
            xlabel('E (\times 10^3 J kg^{-1})')
            ylabel('\zeta')
            grid on
            
            subplot(1,3,2)
            hold on
            plot(enthB_analy_T-273.15, enthB_analy_z, 'r', 'linewidth', 3)
            plot(T(:,1)-273.15, zeta, '-k', 'linewidth', 1.5)
            set(gca, 'xtick', -3:0.5:0.5)
            set(gca, 'ytick', 0:0.1:1)
            set(gca, 'YTickLabel', [])
            
            xlim([-3 0.5])
            set(gca, 'xtick', [-3 -2 -1 0])
            xlabel('T (^\circC)')
            grid on
            
            subplot(1,3,3)
            hold on
            plot(enthB_analy_omega*100, enthB_analy_z, 'r', 'linewidth', 3)
            plot(omega(:,1)*100, zeta, '-k', 'linewidth', 1.5)
            set(gca, 'xtick', 0:0.5:3)
            set(gca, 'ytick', 0:0.1:1)
            set(gca, 'YTickLabel', [])
            
            xlim([-0.5 3])
            xlabel('\omega (%)')
            grid on
            
            set(gcf, 'position', [239 153 1005 412])
        end
        
    end
end