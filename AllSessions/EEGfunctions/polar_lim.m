function h = polar_lim(theta, rho, max_lim, varargin)
% 
%    t = 0:.01:2*pi;
%    polar_lim(t,sin(2*t).*cos(2*t),1.5,'--r')

% max_lim
x_fake=[0 max_lim 0 -max_lim]; 
y_fake=[max_lim 0 -max_lim 0]; 
h_fake=polar(x_fake,y_fake); 
set(h_fake,'Visible','off') 

hold on
h = polar(theta,rho,varargin{:}); 
hold off
end
