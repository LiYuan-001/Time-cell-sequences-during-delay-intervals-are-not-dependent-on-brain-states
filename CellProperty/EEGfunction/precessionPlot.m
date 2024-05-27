function precessionPlot(SpkPhase,SpkPosition,Circ,lin,subPos)

SpkPosition =[SpkPosition';SpkPosition'];
SpkPhase =[SpkPhase';SpkPhase'+360];
gridSize=10;

if isempty(SpkPosition)
    gridEval = zeros(9,9);
    GridEval = zeros(9,9);
    x = 1:10;
    y = 1:10;
else
    
xx=SpkPosition;
yy=SpkPhase;
colormap;
x=linspace(min(xx),max(xx),gridSize);
y=linspace(min(yy),max(yy),gridSize);
gridEval = zeros(length(x)-1,length(y)-1);
for cnt_x=1:length(x)-1
    for cnt_y=1:length(y)-1
        x_ind=intersect(find(xx>x(cnt_x)),find(xx<=x(cnt_x+1)));
        xy_ind=intersect(find(yy(x_ind)>y(cnt_y)), find(yy(x_ind)<=y(cnt_y+1)));
        gridEval(cnt_y, cnt_x)=length(xy_ind);
    end
end
MaxColoum = max(gridEval,[],1);
m=MaxColoum;
for k=1:length(y)-2
    MaxColoum=[MaxColoum;m];
end

MaxColoum(MaxColoum == 0) = 1;
GridEval=gridEval./MaxColoum;
end

% start plotting
subplot(3,5,subPos(1))
plot(SpkPosition,SpkPhase,'.r');hold on;
ylim([0 720]);
xxx=0:1;
%     plot(xxx,180*(2*pi*Wholealpha.*xxx+Wholephi0)/pi);
xlabel('Distance on Arm','Fontsize',10);
ylabel('Theta Phase','Fontsize',10);

CircularMeasure = sprintf('%s%1.2f%s%1.2f%s%1.2f%s%1.2f','Circular Linear fit Alpha=',Circ.Alpha,' Phi0=',180*Circ.Phi0/pi,' Coeff=',Circ.Coeff,' P=',Circ.pValue);
LinearMeasure = sprintf('%s%1.2f%s%1.2f%s%1.2f','Linear fit Alpha=',lin.Alpha,' Phi0=',lin.Phi0,' Coeff=',lin.r);
if Circ.pValue<=0.05
    title({'Theta precession';CircularMeasure;LinearMeasure},'Color','r','Interpreter','None','Fontsize',8)
else
    title({'Theta precession';CircularMeasure;LinearMeasure},'Interpreter','None','Fontsize',8)
end
hold off

subplot(3,5,subPos(2))
surf((x(1:end-1)+ x(2:end))/2,(y(1:end-1)+y(2:end))/2,gridEval); view(2);
shading interp;  hold on;
xlabel('Distance on Arm','Fontsize',10);
ylabel('Theta Phase','Fontsize',10);
title('Spike phase-position density map','Fontsize',8)
h1 = gca;
% add box line
axes(h1);
axis off
hold off

subplot(3,5,subPos(3))
surf((x(1:end-1)+ x(2:end))/2,(y(1:end-1)+y(2:end))/2,GridEval); view(2);
shading interp;  hold on;
xlabel('Distance on Arm','Fontsize',10);
ylabel('Theta Phase','Fontsize',10);
title('Spike phase-position density map-normalized','Fontsize',8)
h2 = gca;
axes(h2);
axis off
hold off
end