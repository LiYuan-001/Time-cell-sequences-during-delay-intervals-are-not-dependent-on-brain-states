function ptempl = parsingtemplate(pmode,mazetype)
%PARSINGTEMPLATE Define zones used to assign zone index to path data
%
% ptempl = PARSINGTEMPLATE(pmode)
%    pmode can be 'fig8:rdscr' (return-delay-stem-choice-reward)
%                 'fig8mouse:rdscr' (return-delay-stem-choice-reward)

% split = strsplit(pmode, ':');
% mazetype = split{1};
% pmode = split{2};
pmode = 'rdscr';
switch mazetype
	case {'fig8', 'fig8rat','fig8Treadmill'}
		if strcmpi(pmode, 'rdscr') % return-delay-stem-choice-reward
			ptempl(1).x = [-70 -70 -25 -25 -70 NaN 25 25 70 70 25];    % return left and right
			ptempl(1).y = [-90 60 60 -90 -90 NaN -90 60 60 -90 -90]; % return left and right
			ptempl(1).zone = 'return';
			
			ptempl(2).x = [-10 -10 10 10 -10]; % delay
			ptempl(2).y = [90 45 45 90 90];   % delay
			ptempl(2).zone = 'delay';
			
			ptempl(3).x = [-20 -20 20 20 -20]; % stem
			ptempl(3).y = [45 -60 -60 45 45]; % stem
			ptempl(3).zone = 'stem';
			
			ptempl(4).x = [-10 -10 10 10 -10]; % choice
			ptempl(4).y = -[30 90 90 30 30];    % choice
			ptempl(4).zone = 'choice';
			
			ptempl(5).x = [-60 -60 -10 -10 -60 NaN 10 10 60 60 10]; % reward
			ptempl(5).y = -[60 90 90 60 60 NaN 60 90 90 60 60];      % reward
			ptempl(5).zone = 'reward';
		end
		
		for i = 1:length(ptempl)
			[ptempl(i).x, ptempl(i).y] = poly2cw(ptempl(i).x, -ptempl(i).y);
		end
		
end

switch mazetype
	case 'fig8mouse'
		for i = 1:length(ptempl)
			[ptempl(i).x, ptempl(i).y] = deal(ptempl(i).x*.5, ptempl(i).y*.5);
		end
end