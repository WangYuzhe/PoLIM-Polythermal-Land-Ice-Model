function u_s = smooth_velocities(u_s, para)
% Apply mild smoothing only to ice areas
% by DeepSeek 2025/06/10

global H_s

ice_mask = H_s > para.Hmin;
for i=find(ice_mask)
    if i>1 && i<length(ice_mask)
        u_s(:,i) = 0.25*u_s(:,i-1) + 0.5*u_s(:,i) + 0.25*u_s(:,i+1);
    end
end
% Ensure zeros in off-ice
u_s(:,~ice_mask) = 0;
end