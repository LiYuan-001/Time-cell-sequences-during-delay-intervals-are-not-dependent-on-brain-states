function sse=myfit(params,Input,Actual_Output)

% Load parameters
A = params(1);
Mu = params(2);
C = params(3);

% sinusoidal fitting curve
Fitted_Curve = A.*cos((Input-Mu)*pi/180)+C; 

% Residuals
Error_Vector = Fitted_Curve - Actual_Output;

% When curvefitting, a typical quantity to
% minimize is the sum of squares error
% You could also write sse as
% sse=Error_Vector(:)'*Error_Vector(:);
sse=sum(Error_Vector.^2);
end
