function mm = rtoff(logicalVec)
%RTON(LOGICALARRAY) (R)e(t)ain off-to-(on) transitions
%
% in a vectors of 0's and 1's retrains the first 1 of each series of
% adjacent 1's
%
% Member of (L)ogical (A)rray (U)titlies (LAU) pack

if isempty(logicalVec)
	mm  = [];
	return;
end

buffer = logicalVec == 1;
buffer(end+1) = 0;
mm = logicalVec & diff(buffer);
% warning('LAU:rtoff', 'rtoff and rton have been substituted on 3-24-2015, make sure the code that''s using them is aware of this!')