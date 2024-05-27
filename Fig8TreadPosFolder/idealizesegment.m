function [x, y] = idealizesegment(x0,y0,xNode,yNode,entryExit)

% comments added by @author: Sia @date 9/15/15 5:27 PM
% modified code to avoid the new behavior of @minmax in Matlab 2016b
% 4/6/2017 (used to be  mean(minmax(yNode)) as well as for xNode)

if all(~ismember({'l','r'},entryExit)) % a vertical arm @author: Sia @date 9/15/15
    x = x0;
    y = y0;        
   % straightens the path vertically by replacing the horizontal movement by a constant
    if nanmean(x) < -30
        x(:) = mean([min(xNode(:)), max(xNode(:))])+diff([min(xNode(:)), max(xNode(:))])/3;
    elseif nanmean(x) > 30
        x(:) = mean([min(xNode(:)), max(xNode(:))])-diff([min(xNode(:)), max(xNode(:))])/3;
    else
        x(:) = mean([min(xNode(:)), max(xNode(:))]);
    end

elseif all(~ismember({'t','b'},entryExit)) % a base arm @author: Sia @date 9/15/15
    x = x0;
    y = y0;
    
    if nanmean(y) < -50
        y(:) = mean([min(yNode(:)), max(yNode(:))])+ diff([min(yNode(:)), max(yNode(:))])/3;
    elseif nanmean(y) > 50
        y(:) = mean([min(yNode(:)), max(yNode(:))])- diff([min(yNode(:)), max(yNode(:))])/3;
    else
        y(:) = mean([min(yNode(:)), max(yNode(:))]); % gets rid of vertical displacement
    end
else
    if ismember('t',entryExit)
        dy = max(yNode);
    else
        dy = min(yNode);
    end
    if ismember('l',entryExit)
        dx = min(xNode);
    else
        dx = max(xNode);
	end
	% project path onto a quarter circle:
    a = cart2pol(x0-dx,y0-dy);
    r = range(xNode)/2; % half length node square
    [x, y] = pol2cart(a,r);
    x = x+dx;
    y = y+dy;
end
end