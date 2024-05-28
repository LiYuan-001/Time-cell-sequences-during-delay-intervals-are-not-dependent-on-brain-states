function PlotThetaHist(n,Phase,spikeRadPhase,xDetail,p,MinH,MaxH,pVal,MaxSpike,x)
if pVal<0.05
    bar([x x+360],[n n],1,'r')
    title({'p<0.05';'von Mises fitting'})
else
    bar([x x+360],[n n],1,'b')
    title('von Mises fitting');
end
% Curve fitting
if p.fit && ~isempty(Phase)
    % von Mises fitting of original data
    [mu kappa] = circ_vmpar(spikeRadPhase);
    % get von Mises distribution
    vm = circ_vmpdf(circ_ang2rad(xDetail),mu,kappa);
    vm = vm/sum(vm)*length(xDetail)/(360/p.degBin);
    vm=vm(:);
    %plot
    hold on
    plot([xDetail xDetail+360],[vm' vm'],'k-')
    %     [A,Mu,C] = sinusoidalFitting(nT,x);
    %     y = fittedCurve(xDetail,A,Mu,C);
    %     hold on
    %     plot([xDetail xDetail+360], [y y],'k-')
end
% Spike raster
if p.raster && ~isempty(Phase)
    idx = min(length(Phase), MaxSpike);
    plot([Phase(1:idx) Phase(1:idx)],[MaxH-0.01 MaxH],'r')
    plot([Phase(1:idx)+360 Phase(1:idx)+360],[MaxH-0.01 MaxH],'r')
end
% Setting
set(gca,'XTick',[0 180 360 540 720]); axis([0 720 MinH MaxH]);
set(gca,'Box','Off')
set(gca, 'TickDir', 'Out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca,'FontSize',8);
set(gca,'XTick',[0 180 360 540 720])

end
