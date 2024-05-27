function mm = rton(logicalVec)
%RTOFF(LOGICALARRAY) (R)e(t)ain on-to-(off) transitions
%
% in a vectors of 0's and 1's retains the last 1 of each series of
% adjacent 1's
%
% Member of (L)ogical (A)rray (U)titlies (LAU) pack

mm = logicalVec(:)'==1 & diff([0 logicalVec(:)']==1);
% warning('LAU:rton', 'rtoff and rton have been substituted on 3-24-2015, make sure the code that''s using them is aware of this!')

if iscolumn(logicalVec)
	mm = mm';
end