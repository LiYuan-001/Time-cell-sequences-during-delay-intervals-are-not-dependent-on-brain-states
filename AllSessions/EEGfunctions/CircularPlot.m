function CircularPlot(Num,spikeRadPhase,r,ang,p,Max)
if Num ~= 0
    % polar_lim(spikeRadPhaseT,ones(length(ts),1)*MaxT,'o')
    [tout,rout]=rose(spikeRadPhase,360/p.degBin);
    polar_lim(tout,rout./Num,Max);
    hold on
end
h=compass_lim(r*cos(circ_ang2rad(ang)),r*sin(circ_ang2rad(ang)),Max);
set(h,'Color',[1 0 0],'LineWidth',4)
ylabel('Spike phase')
Xlabel = sprintf('%s%u','Spike number: ',Num);
xlabel(Xlabel);
end
