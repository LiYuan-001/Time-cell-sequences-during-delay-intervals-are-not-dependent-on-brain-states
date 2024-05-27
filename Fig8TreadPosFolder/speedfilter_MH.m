function [x,y] = speedfilter_MH(x,y,t,p)
if ~isfield(p,'pix2cm')
    pix2cm = 0.5;
else
    pix2cm = p.pix2cm;
end
a = hypot(pix2cm*diff(x),pix2cm*diff(y));
a = [a(1),a];
pertime = diff(t)*10^-6;  %conversion to seconds 
pertime = [pertime(1),pertime];
rates = a./pertime;
filter_values = find(rates>p.speed);
%Use this to plot
figure
scatter(x(filter_values),y(filter_values),'r*');hold on
plot(x,y,'b'); hold off
title('Before smooth & speed filtering');
% Filter now points marked in red
x(filter_values) = NaN;
y(filter_values) = NaN;
end