function [ansr] = CCCGraphEval(dt, Ref_rho, alfa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is the version 2 of code
%   Compute the Concordance Correlation Coefficient proposed by Dr. Lin,
%   L.K. 
%   Reference: Lin, L.K., 1989, A Concordence Correlation Coefficient to Evaluate
%   Reproducibility, Biometrics, 45:255-268
%   Confidence Interval is corrected on 3/11 by his correction paper in
%   2000, Biometrics 56 (1) 324-325

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
else
end

%   1. Point estimation of CCC
[nn mm] = size(dt);         ansr    = CCC(dt, alfa) ;

%   2. Graph in (Mean, Difference)
%   2.1 B-A method
DiffX       = dt(:,2) - dt(:,1);
dt_t        = [mean(dt')' DiffX];

d_bar       = mean(DiffX);
BA_limit    = tinv(1-alfa/2, nn-1) * std(DiffX);
BA_upr      = d_bar + BA_limit;      
BA_lwr      = d_bar - BA_limit;

ansr.d_bar  = d_bar;                ansr.BA_limit   = BA_limit;
ansr.Sd     = std(DiffX) ;

%   2.2 Reference band
ansr.S  = ansr.Sd / sqrt(2*(1-ansr.rho)) ;
Wr      = tinv(1-alfa/2, nn-1) * ansr.S * sqrt( 2*(1-Ref_rho) ) ;

Xmax    = floor(max(dt_t(:,1)))+1.5 ;          
Xmin    = floor(min(dt_t(:,1)))-0.5 ;
Ymin    = floor(min(min(dt_t(:,2)), min(-Wr, -BA_limit))) ;
Ymax    = floor(max(max(dt_t(:,2)), max(Wr, BA_limit))) + 1 ;
YYmax   = max(abs(Ymin), abs(Ymax)) ;

jrc_1   = dt_t(:,2) >= -Wr;             % Judge by Reference band
jrc_2   = dt_t(:,2) <=  Wr;
jr_ref  = jrc_1.*jrc_2;              
nrc     = sum(jr_ref);                  % # of data within reference band

ansr.Wrb        = Wr ;
ansr.WithinBand = nrc ;
ansr.Outliers   = nn - nrc ;
ansr.N          = nn ;

S1 = std(dt(:,1))^2 ;               S2 = std(dt(:,2))^2 ;
p1 = 2*fcdf(S1/S2, nn-1, nn-1);     p2 = 2*fcdf(S1/S2, nn-1, nn-1, 'upper') ;
ansr.F_test   = min([p1 p2]) ;

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

lgd = legend(['RB with CCC = ' num2str(Ref_rho,  '%3.2f')], 'Limits of Agreement');
set(lgd,'EdgeColor',[1 1 1], 'FontSize', 10);

text(1.1*Xmin, -0.9*YYmax, ['Estimates: CCC = ' num2str(ansr.point,  '%5.3f') ', \rho = ', num2str(ansr.rho, '%5.3f')...
    ', C_b = ' num2str(ansr.Cb, '%5.3f') ', and # of outliers = ', num2str(ansr.Outliers, '%3.0f')],...
        'HorizontalAlignment','left', 'Fontsize', 9);
%     
% if ansr.F_test < 0.05
%     text(1.1*Xmin, -0.7*YYmax, {'Reference Band may not be valid since two variances are signficantly different'}, ...
%         'Fontsize', 10, 'Color', 'red');
% else
% end    
    
axis([Xmin Xmax -YYmax/0.9 YYmax/0.9]);
xlabel('Mean, ($X_1$ + $X_2$)/2', 'Interpreter','latex', 'fontsize', 11) ;
% ylabel('$\sqrt{n}$ ($X_2$ - $X_1$)', 'Interpreter','latex', 'fontsize', 11) ;
ylabel('Difference ($X_2$ - $X_1$)', 'Interpreter','latex', 'fontsize', 11) ;
box off
hold off
%   Final revision on February 16, 2021 by Jongphil Kim