function mu = circ_mean2(DataOri, tOri)
% Computes the mean direction for circular data given as a distribution.
%     Data  value at each angle
%     t		angular bin (ex. 15:30:360)
%
% Takuma Kitanishi, CBM, NTNU, 2011

% Column vectorize
Data = DataOri(:);
t = tOri(:);

% ang2rad
t = t * pi /180;

% compute weighted sum of cos and sin of angles
r = Data'*exp(1i*t);

% obtain mean by
mu = angle(r);

% rad2angle
mu = mu / pi *180;

% 0 to 360
mu = mod(mu,360);
end
