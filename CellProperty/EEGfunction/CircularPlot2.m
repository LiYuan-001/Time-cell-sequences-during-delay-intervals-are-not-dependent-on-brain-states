function CircularPlot2(Num,spikeRadPhase,r,ang,p,pVal,Max)
if Num ~= 0
    % polar_lim(spikeRadPhaseT,ones(length(ts),1)*MaxT,'o')
    [tout,rout]=rose(spikeRadPhase,360/p.degBin);
    polar_lim(tout,rout./Num,Max);
    hold on
end
h=compass_lim(r*cos(circ_ang2rad(ang)),r*sin(circ_ang2rad(ang)),Max);
if pVal <= 0.05
    set(h,'Color',[1 0 0],'LineWidth',3)
else
    set(h,'Color',[0 0 1],'LineWidth',3)
end
ylabel('Spike phase')
Xlabel = sprintf('%s%u','Spike number: ',Num);
xlabel(Xlabel);
end
