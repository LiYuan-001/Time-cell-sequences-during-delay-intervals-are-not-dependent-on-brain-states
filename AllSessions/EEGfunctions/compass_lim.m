function h = compass_lim(x, y, max_lim)
% 
%   x = eig(randn(20,20));
%   compass_lim(x,10);
%   %% or
%   compass_lim(real(x),imag(x),10);

if nargin == 2
    max_lim = y;
    y = imag(x);
    x = real(x);
end

% max_lim ‚ÌƒTƒCƒY‚Å•\¦
x_fake=[0 max_lim 0 -max_lim]; 
y_fake=[max_lim 0 -max_lim 0]; 
h_fake=compass(x_fake,y_fake); 
set(h_fake,'Visible','off') 

% ’Ê?í‚Ìƒv?ƒbƒg‚ğ?d‚Ë?‘‚?
hold on
h = compass(x,y); 
hold off
end
