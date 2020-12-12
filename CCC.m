function [Concordance Accuracy Precision] = CCC(datas, alfa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Concordance Correlation Coefficient proposed by Dr. Lin, L.K. 
%   Reference: Lin, L.K., 1989, A Concordence Correlation Coefficient to Evaluate
%   Reproducibility, Biometrics, 45:255-268
%
%   Confidence Interval is corrected on 3/11 by his correction paper in
%   2000, Biometrics 56 (1) 324-325
%  
%   coded on February 19, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Default Settings and Input variables

if nargin == 1          alfa = 0.05;    else end
type_x = 1;

% xx = datas.dt(:,1);         yy = datas.dt(:,2);
% VarName = datas.name;
% if isempty(VarName)
%     VarName = {'Mean(X_1 + X_2) ', 'Difference (X_2 - X_1)'};
% else
% end
% type_x = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Point Estimate

xx = datas(:,1);        yy = datas(:,2);
Non_nan = ~ (isnan(xx) | isnan(yy));    %   Check if there are any missing values
dt = [xx(Non_nan) yy(Non_nan)];
[nn mm] = size(dt);

MeanX = mean(dt(:,1));      MeanY = mean(dt(:,2));
VarCov = cov(dt, 1);                    %   Scaled by n, not n-1
VarX = VarCov(1,1);         VarY = VarCov(2,2);
Sxy = VarCov(1,2);

rr = 2*Sxy / (VarX + VarY + (MeanY - MeanX)^2);     %   Concordance
rho_s       = corrcoef(dt(:,1), dt(:,2));     
Pearson_r   = rho_s(1,2);                           % Pearson Correlation Coefficient
Cb          = rr / Pearson_r ;

nu = (MeanX - MeanY)/sqrt(sqrt(VarX)*sqrt(VarY));
omega = sqrt(VarY/VarX);
Ca = 2/(omega + 1/omega + nu^2);                    %   Accuracy

r = Sxy/sqrt(VarX*VarY);                            %   Precision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Confidence Limits

Zrr = 0.5*log((1 + rr)/(1 - rr));                   %   Transformation
Lca = log ( Ca/(1 - Ca) );
Zr = 0.5*log((1 + r)/(1 - r));

if type_x == 0              %   when X is random    
    %   1) Concordance
    VarZ1 = (1 - r^2) * rr^2 / ((1 - rr^2) * r^2);
    VarZ2 = 2 * rr^3 * (1 - rr) * nu^2 / (r * (1 - rr^2)^2);
    VarZ3 = rr^4 * nu^4 / (2 * r^2 * (1 - rr^2)^2);
    VarZrr = (VarZ1 + VarZ2 - VarZ3)/(nn - 2);  
    
    %   2) Accuracy
    VarLca1 = Ca^2 * nu^2 * (omega + 1/omega - 2*r) + 0.5 * Ca^2 * ...
        (omega^2 + 1/omega^2 + 2 * r^2) + (1 + r^2)*(Ca * nu^2 - 1);
    VarLca2 = (nn - 2) * (1 - Ca)^2;
    VarLca = VarLca1 / VarLca2;
    
    %   3) Precision
    VarZr = 1/(nn - 3);

elseif type_x == 1          %   when X is fixed
    %   1) Concordance
    VarZ1 = (1 - r^2) * rr^2 / ( (nn-2) * (1-rr^2)^2 * r^2 );
    VarZ2 = omega * nu^2 * rr^2 + (1 - rr * r * omega)^2 + ...
        0.5 * omega^2 * rr^2 * (1 - r^2);
    VarZrr = VarZ1 * VarZ2;

    %   2) Accuracy
    VarLca1 = nu^2 * omega * Ca^2 * (1 - r^2) + 0.5 * (1 - omega*Ca)^2 *...
        (1 - r^4);
    VarLca2 = (nn - 2) * (1 - Ca)^2;
    VarLca = VarLca1 / VarLca2;
    
    %   3) Precision
    VarZr = (1 - 0.5 * r^2)/(nn - 3);
else
    error('Type of Xs should be 0 or 1');
end

%   Concordance
Lwr1Zrr = Zrr - norminv(1-alfa) * sqrt(VarZrr);         %   1-sided
Lwr2Zrr = Zrr - norminv(1-alfa/2) * sqrt(VarZrr);       %   2-sided
Upr2Zrr = Zrr + norminv(1-alfa/2) * sqrt(VarZrr);
    
Lwr1rr = (exp(2*Lwr1Zrr) - 1)/(exp(2*Lwr1Zrr) + 1);
Lwr2rr = (exp(2*Lwr2Zrr) - 1)/(exp(2*Lwr2Zrr) + 1);
Upr2rr = (exp(2*Upr2Zrr) - 1)/(exp(2*Upr2Zrr) + 1);

Concordance.point   = rr;
Concordance.lwr1    = Lwr1rr;   %   1-sided lower
Concordance.lwr2    = Lwr2rr;   %   2-sided lower
Concordance.upr2    = Upr2rr;   %   2-sided upper
Concordance.Cb      = Cb ;
Concordance.rho     = Pearson_r ;
Concordance.alfa    = alfa ;

%   Accuracy
Lwr1Lca = Lca - norminv(1-alfa) * sqrt(VarLca);         %   1-sided
Lwr2Lca = Lca - norminv(1-alfa/2) * sqrt(VarLca);       %   2-sided
Upr2Lca = Lca + norminv(1-alfa/2) * sqrt(VarLca);
       
Lwr1Ca = exp(Lwr1Lca)/(1 + exp(Lwr1Lca));
Lwr2Ca = exp(Lwr2Lca)/(1 + exp(Lwr2Lca));
Upr2Ca = exp(Upr2Lca)/(1 + exp(Upr2Lca));

Accuracy.point = Ca;
Accuracy.lwr1 = Lwr1Ca;
Accuracy.lwr2 = Lwr2Ca;
Accuracy.upr2 = Upr2Ca;

%   Precision
Lwr1Zr = Zr - norminv(1-alfa) * sqrt(VarZr);            %   1-sided
Lwr2Zr = Zr - norminv(1-alfa/2) * sqrt(VarZr);          %   2-sided
Upr2Zr = Zr + norminv(1-alfa/2) * sqrt(VarZr);
    
Lwr1r = (exp(2*Lwr1Zr) - 1)/(exp(2*Lwr1Zr) + 1);
Lwr2r = (exp(2*Lwr2Zr) - 1)/(exp(2*Lwr2Zr) + 1);
Upr2r = (exp(2*Upr2Zr) - 1)/(exp(2*Upr2Zr) + 1);

Precision.point = r;
Precision.lwr1 = Lwr1r;
Precision.lwr2 = Lwr2r;
Precision.upr2 = Upr2r;

%   End of Code