function [map pospdf] = ratemap(spkx,spky,posx,posy,post,h,mapAxis1,mapAxis2)
%I'm just removing nans here, but NaNs come from dividing linear track data
%by direction, and hence mark discontinuities in the data
%over which edges really shouldn't be smoothed.  This needs to be taken
%care of properly! But I'm hoping this will do for now.
naninds = isnan(posx+posy+post);
posx(naninds)=[];posy(naninds)=[];post(naninds) = [];

invh = 1/h;
map = zeros(length(mapAxis2),length(mapAxis1));
pospdf = zeros(length(mapAxis2),length(mapAxis1));

yy = 0;
for y = mapAxis2
    yy = yy + 1;
    xx = 0;
    for x = mapAxis1
        xx = xx + 1;
        [map(yy,xx),pospdf(yy,xx)] = rate_estimator(spkx,spky,x,y,invh,posx,posy,post);
    end
end

% Normalize the pospdf (Should actually normalize the integral of the
% pospdf by dividing by total area too, but the functions making use of the 
% pospdf assume that this is not done (see mapstat()).
pospdf = pospdf ./ sum(sum(pospdf)); % Position Probability Density Function
end