function precessionPlotOrigin(SpkPhase,SpkPosition,fieldPos,P,s,b)
if size(SpkPhase,1) == 1
    SpkPhase = [SpkPhase,SpkPhase+360];
    SpkPosition = [SpkPosition,SpkPosition];
else
    SpkPhase = [SpkPhase;SpkPhase+360];
    SpkPosition = [SpkPosition;SpkPosition];
end
% generate coordinates for plotting fitting line
x1 = 0;
x2 = 1;
y1 = 2*pi*x1*s*180/pi+b*180/pi;
y2 = 2*pi*x2*s*180/pi+b*180/pi;

% start plotting
if P <= 0.05
    plot(SpkPosition,SpkPhase,'.r','MarkerSize',2);hold on;
    for i = -3:3
        plot([min(fieldPos),max(fieldPos)],[y1,y2]+i*360,'r')
    end
else
    plot(SpkPosition,SpkPhase,'.b','MarkerSize',2);hold on;
    for i = -3:3
        plot([min(fieldPos),max(fieldPos)],[y1,y2]+i*360,'b')
    end
end
% ylabel('Theta Phase','Fontsize',10);
ylim([0 720])

end