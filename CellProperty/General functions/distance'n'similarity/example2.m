clc, clear, close all

%% generate test data
x = [0 5];
y = [5 0];
                                                              
%% calculate the metrics
[ed, es] = eucmetric(x, y); 
[cd, cs] = cosmetric(x, y);  
[ad, as] = angmetric(x, y);                                    
                          
%% plot the results
figure(1)
plot([0 x(1)], [0 x(2)], '-ob', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot([0 y(1)], [0 y(2)], '--xr', 'MarkerSize', 10, 'LineWidth', 2)
grid minor
xlim([-10 10])
ylim([-10 10])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('\itX')
ylabel('\itY')
legend('vector 1', 'vector 2')

%% show the results
commandwindow
table(categorical({'Euclidian distance'; 'Euclidian similarity'; ...
                   'Cosine distance'; 'Cosine similarity'; ...
                   'Angular distance'; 'Angular similarity'}), ...
                   [ed; es; cd; cs; ad; as], ...
                   'VariableNames', {'Property', 'Value'})