function plotSimulationResult(T,Y,FLUX,xrange,yrange, sampling_time)

global strain_no
[Yidx, Fidx] = setIndex();

[ T_exp, X_exp, GLCex_exp, ACEex_exp ] = ExpDataForBatchCulture();

%% Time course 
scrsz = get(0,'ScreenSize');
figure('Position',[20 scrsz(4)*0.3 scrsz(3)*0.3 scrsz(4)*0.6]);

subplot(2,1,1)
plot(T_exp,[X_exp GLCex_exp*1e-3*180 ACEex_exp*1e-3*60],'o');hold on
plot(T,[ Y(:,Yidx.X) Y(:,Yidx.GLCex)*1e-3*180 Y(:,Yidx.ACEex)*1e-3*60 ],'LineWidth',2);hold on
xlabel('Time (h)','FontSize',10,'FontName','Arial');
ylabel('Concentration (g/L)','FontSize',10,'FontName','Arial');
legend('X','GLCex','ACEex','Location','northeast');
xlim(xrange);
ylim(yrange);
set(gca,'XTick',xrange(1):5:xrange(2));
set(gca,'YTick',0:2:20);
set(gca,'FontSize',10,'FontName','Arial');

switch strain_no
    case 1
        title('Wild Type');
    case 26
        title('\Delta\itpykA/pykF');
    case 4
        title('\Delta\itpgi');
    case 25
        title('\Delta\itppc');
    otherwise
        fprintf('Unexpected Strain Number!\n');
end

%% Comparison between the experimental and simulated fluxes
[Flux_exp, Flux_sim] = rearrange_exp_sim_flux(T, FLUX, sampling_time);
subplot(2,1,2)
x = 1:1:size(Flux_sim,2);
bar(x, [Flux_exp(1,:)' Flux_sim(1,:)'],'grouped');
ax = get(gca);
cat = ax.Children;
xlabel('Reaction index (-)','fontname','arial','fontweight','bold','fontsize',12);
ylabel('Relative flux (-)','fontname','arial','fontsize',12);
ylim([0 Inf]);
xlim([0 29]); xtick=1:1:28; 
title(' 5h','fontname','arial','fontsize',12)
set(cat(2),'FaceColor','g','BarWidth',2);
set(cat(1),'FaceColor','r','BarWidth',2);
set(gcf,'color','white');
legend('Experiment','Simulation');

return
