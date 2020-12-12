function [ansr] = CCCGraphEval(dt, Ref_rho, alfa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is developed to conduct the graphical evaluation based on the CCC, 
%   Concordance Correlation Coefficient proposed by Dr. Lin, L.K. 
%
%   1st column of dts is 1st measurement, 2nd column is 2nd measurement

%   Call CCC.m file to compute CCC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Data arrangement
if nargin == 1      
    Ref_rho = 0.75 ;    alfa = 0.05;
elseif nargin == 2
    alfa = 0.05;
else
end

if sum(sum(isnan(dt))) > 0
    xx = datas(:,1);        yy = datas(:,2);
    Non_nan = ~ (isnan(dt(:,1)) | isnan(dt(:,2)));    %   Check if there are any missing values
    dt = dt(Non_nan,:);
%     disp('Number of NaN in each column');
%     disp(sum(isnan(dt)));
%     error('Please remove NaN in dataset');
else
end

%   1. Point estimation of CCC
[nn mm] = size(dt);                 nu = 2*(nn-1);
CritVal = tinv(1-alfa/2, nu);
ansr = CCC(dt, alfa) ;

%   2. Graph in (Mean, Difference)
%   2.1 B-A method
DiffX = dt(:,2) - dt(:,1);
dt_t = [mean(dt')' DiffX];

d_bar = mean(DiffX);
BA_limit = tinv(1-alfa/2, nn-1) * std(DiffX);
BA_upr = d_bar + BA_limit;      BA_lwr = d_bar - BA_limit;

ansr.d_bar      = d_bar;        ansr.BA_limit   = BA_limit;
ansr.Sd         = std(DiffX) ;

%   2.2 Reference band
S       = sqrt( nn*(std(dt(:,1),1)^2 + std(dt(:,2),1)^2)/(2*(nn-1)) ) ;
Wr      = S * CritVal * sqrt(2*(1-Ref_rho));

Xmax    = floor(max(dt_t(:,1)))+1.5 ;          
Xmin    = floor(min(dt_t(:,1)))-0.5 ;
Ymin    = floor(min(min(dt_t(:,2)), min(-Wr, -BA_limit))) ;
Ymax    = floor(max(max(dt_t(:,2)), max(Wr, BA_limit))) + 1 ;
YYmax   = max(abs(Ymin), abs(Ymax)) ;

jrc_1   = dt_t(:,2) >= -Wr;          % Judge by Reference band
jrc_2   = dt_t(:,2) <=  Wr;
jr_ref  = jrc_1.*jrc_2;              
nrc     = sum(jr_ref);                  % # of data within reference band

ansr.S          = S;
ansr.Wr         = Wr ;
ansr.WithinBand = nrc ;
ansr.Outliers   = nn - nrc ;
ansr.N          = nn ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2 ;
figure
p1 = scatter(dt_t(:,1), dt_t(:,2), 'r', 'filled');
p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on

p2= plot([Xmin Xmax], [0 0], 'linewidth', 1.5, 'Color', 'k');
p2.Annotation.LegendInformation.IconDisplayStyle = 'off';

p3 = plot([Xmin Xmax], -Wr*ones(1,2) , '-',  'linewidth', lw, ...
    'Color', [0.494117647058824 0.184313725490196 0.556862745098039]) ;
p4 = plot([Xmin Xmax], Wr*ones(1,2) , '-',  'linewidth', lw, ...
    'Color', [0.494117647058824 0.184313725490196 0.556862745098039]);
p4.Annotation.LegendInformation.IconDisplayStyle = 'off';

p5 = plot([Xmin Xmax], BA_upr*ones(1,2) , 'b-.',  'linewidth', lw);
p6 = plot([Xmin Xmax], BA_lwr*ones(1,2) , 'b-.',  'linewidth', lw);
p6.Annotation.LegendInformation.IconDisplayStyle = 'off';

lgd = legend(['Reference band with CCC = ' num2str(Ref_rho,  '%3.2f')], 'Limit of Agreement by Bland-Altman');
set(lgd,'EdgeColor',[1 1 1], 'FontSize', 10);

text(1.1*Xmin, -0.9*YYmax, ['Estimates: CCC = ' num2str(ansr.point,  '%5.3f') ', \rho = ', num2str(ansr.rho, '%5.3f')...
    ', C_b = ' num2str(ansr.Cb, '%5.3f') ', and # of outliers = ', num2str(ansr.Outliers, '%3.0f')],...
        'HorizontalAlignment','left', 'Fontsize', 9);
    
axis([Xmin Xmax -YYmax/0.9 YYmax/0.9]);
xlabel('Mean, (X_1 + X_2)/2', 'fontsize',11);      
ylabel('Difference (X_2 - X_1)', 'fontsize',11);
box off
hold off
%   Final revision on February 20, 2020 by Jongphil Kim