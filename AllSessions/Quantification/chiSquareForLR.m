% chi-square test foe left right significant cell and assemblies
sessionName = {'on10','off10','on30','off30'};
armName = {'return','delay','stem','choice'};

cellNum = [176,193,147,148];
assemblyNum = [99,83,69,77];

cellMat = [98,10,12,77;97,9,18,88;99,11,14,81;95,11,18,86];
assemblyMat = [42,8,5,20;49,5,5,22;49,9,4,24;37,4,7,19];
    
for j = 1:length(armName)
    for k = 1:length(sessionName)-1
        for m = k+1:length(sessionName)
            testName = sprintf('Cell_%s_%s_%s',armName{j},sessionName{k},sessionName{m});
            stat.(testName) = chi_Independence(cellMat(k,j),cellNum(j),cellMat(m,j),cellNum(j));
            
            testName = sprintf('Assembly_%s_%s_%s',armName{j},sessionName{k},sessionName{m});
            stat.(testName) = chi_Independence(assemblyMat(k,j),assemblyNum(j),assemblyMat(m,j),assemblyNum(j));
        end
    end    
end

for k = 1:length(sessionName)
    for j = 1:length(armName)
        testName = sprintf('Cell_Assembly_%s_%s',armName{j},sessionName{k});
        stat.(testName) = chi_Independence(cellMat(k,j),cellNum(j),assemblyMat(k,j),assemblyNum(j));
    end
end

stat

function stat = chi_Independence(n1,N1,n2,N2)
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2);
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
stat.chi2stat = sum((observed-expected).^2 ./ expected);
stat.p = 1 - chi2cdf(stat.chi2stat,1);
stat.Input = [n1,N1,n2,N2];
end