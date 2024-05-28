function Deriv_2 = Deriv(waveForm)

[~,troughID] = min(waveForm);

if length(waveForm)-8 > troughID
    curve = waveForm(troughID+3:end-5);
else
    curve = waveForm(troughID+3:end);
end

a=diff(curve);
b=diff(a);
Deriv_2 = mean(b);
end