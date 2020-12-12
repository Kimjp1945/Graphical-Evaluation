con%%  Github version
clear all;
ba = xlsread('real_example_data.xlsx', 'BA');
[Manual, header1]   = xlsread('real_example_data.xlsx', 'manual');
[Ensemble, header2] = xlsread('real_example_data.xlsx', 'ensemble');


[a] = CCCGraphEval([ba(:,4) ba(:,2)]);
axis([215 650 -200 200]);
ylabel('Difference, (Large meter - Mini meter)');
title('PEFR Data', 'fontsize', 12);

TestOrder = [ones(1,32); 2*ones(1,32)];
TestOrder = TestOrder(:);
Manual1     = Manual(TestOrder == 1, :) ;
Manual2     = Manual(TestOrder == 2, :) ;
Ensemble1   = Ensemble(TestOrder == 1, :);
Ensemble2   = Ensemble(TestOrder == 2, :);

%%%%%%%%%%%%%%%
ShortLong = 1 ;
[a3] = CCCGraphEval(log([Manual1(:, ShortLong) Manual2(:, ShortLong)])) ;
title(['Manual, ' header1{1,ShortLong+1}]);  

[a4] = CCCGraphEval(log([Ensemble1(:, ShortLong) Ensemble2(:, ShortLong)])) ;
title(['Ensemble, ' header1{1,ShortLong+1}]); 

%%%%%%%%%%%%%%%
Vol = 2 ;
[a3] = CCCGraphEval(log([Manual1(:, Vol) Manual2(:, Vol)])) ;
title(['Manual, ' header1{1,Vol+1}]);  

[a4] = CCCGraphEval(log([Ensemble1(:, Vol) Ensemble2(:, Vol)])) ;
title(['Ensemble, ' header1{1,Vol+1}]); 

%   End of code