function possibleArms = armsbetween(arm1, arm2)
allArms = {'A16';
	'A12';
	'A25';
	'A23';
	'A45';
	'A34'};
arm1nodes = finite(str2double(num2cell(arm1)));
arm2nodes = finite(str2double(num2cell(arm2)));
possibleArms = [arm1nodes(1), arm2nodes(1);
	arm1nodes(1), arm2nodes(2);
	arm1nodes(2), arm2nodes(1);
	arm1nodes(2), arm2nodes(2);];
possibleArms = [min(possibleArms, [], 2), max(possibleArms, [], 2)];
possibleArms = strcat('A', num2str(possibleArms, '%d%d'));
possibleArms = row2cell(possibleArms);
possibleArms = intersect(possibleArms, allArms);
