function [Phase, Amp] = thetaPhase2(Eeg)

% Hilbert transform
Z = hilbert(Eeg);

% Wave amplitude
Amp = abs(Z);

% EEG phase in rad
Phase = angle(Z);

% Rad to Degree (-180 to +180)
Phase = Phase / pi *180;

% Degree (0 to +360)
Phase = Phase + 180;
end