clc, clear, close all

%% generate test data
x = [1 0 1 0;
     2 5 2 5];
y = [4 5 -4 5;
     8 0 -8 5];
 
%% calculate the metrics
[ed, es] = eucmetric(x, y); 
[cd, cs] = cosmetric(x, y);  
[ad, as] = angmetric(x, y);  

%% show the results
commandwindow 
table(x, y, 'VariableNames', {'x', 'y'})
table(categorical({'Euclidian distance'; 'Euclidian similarity'; ...
                   'Cosine distance'; 'Cosine similarity'; ...
                   'Angular distance'; 'Angular similarity'}), ...
                   [ed; es; cd; cs; ad; as], ...
                   'VariableNames', {'Property', 'Value'})
