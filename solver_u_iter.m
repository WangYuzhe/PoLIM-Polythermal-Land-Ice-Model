function [u_s, u, visc_s, visc, strainHeat] = solver_u_iter(visc_s,visc,AGlen_s,para)
% create on 2023/12/1
% based on De Smedt et al., 2010, JG

global iter_u Ms u_s_lst

iter_u = 0;
while 1
    iter_u = iter_u + 1;
    fprintf('iter_u: %d\n', iter_u)

    [u_s] = solver_u_core(visc_s, visc, AGlen_s, para);
    %----------------------Begin Picard iteration----------------------
    if iter_u>2
        u_s_now = u_s;
        Cs = u_s_now - u_s_lst;
        Sita = acos(Cs'*C/(sumsqr(Cs)*sumsqr(C)));
        if isequal(Sita<=pi/8, ones(Ms,Ms))
            mu1 = 2.5;
        elseif isequal(Sita>pi/8,ones(Ms,Ms)) && isequal(Sita<19*pi/20, ones(Ms,Ms))
            mu1 = 1;
        elseif isequal(Sita>=19*pi/20, ones(Ms,Ms))
            mu1 = 0.5;
        end
        u_s = u_s_lst + mu1*Cs;

        if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now)<1e-4
            break
        end
    end

    if iter_u>1
        C = u_s - u_s_lst;
    end

    u_s_lst = u_s;

    if iter_u>=para.iter_max
        break
    end
    %-----------------------End Picard iteration-----------------------
    u = staggerX2main(u_s);
    [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, para);
end

end