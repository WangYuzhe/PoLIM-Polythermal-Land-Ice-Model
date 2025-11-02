% plot transient results as movie

xx = ones(N,1)*xi/1e3;

zmax = ceil(max(hS)) + 100;
zmin = floor(min(hB)) - 100;

for i=1:numTimeStep
    figure(1), clf(1)
    yy = ones(N,1)*hB + zeta*At_H(i,:);
    ui = At_u(:,:,i);
    umax = max(max(ui));
    umin = 0;
    uint = (umax - umin)/10;

    [max_row, max_col] = wyz_max_index(ui);
    [min_row, min_col] = wyz_min_index(ui);
    
    hold on
    contourf(xx,yy,ui,30,'LineStyle', 'none')
    plot(xi(max_col)/1e3, hB(max_col) + (max_row-1)/(N-1)*H(max_col), 'ro');
    plot(xi(min_col)/1e3, hB(min_col) + (min_row-1)/(N-1)*H(min_col), 'b^');
    
    ch = colorbar('EastOutside');
    set(get(ch,'title'),'string','[m/a]','FontSize',14);
    colormap('Jet')
    
    plot(xi/1e3,hB,'k','linewidth',1.5)
    plot(xi/1e3,At_hS(i,:),'k','linewidth',1.5)
    hold off
    
    title(['Model Year: ',num2str(arrayTime(i),'%10.0f')],'FontSize',16)
    
    ylim([zmin, zmax])
    
    xlabel('x (km)','Fontsize',14)
    ylabel('z (m)','Fontsize',14)
    set(gca,'FontSize',14)
    set(gcf,'position',[200,200,800,400])
    box on
    
    F(i) = getframe(gcf);
end

outfile = './outputs/fixed30_ideal_glac.avi';
videoFile = VideoWriter(outfile);
videoFile.FrameRate = 4;
open(videoFile);
writeVideo(videoFile,F);
close(videoFile);

