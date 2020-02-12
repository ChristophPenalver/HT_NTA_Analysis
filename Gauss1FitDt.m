function [fitresult,gof,locations] = Gauss1FitDt(edges, histdata,output,Name)
%CREATEFIT(EDGES,H4SMOOTHED)
%  Create a fit.
%% Fit: 'Gaus2FitofSize'.
[xData, yData] = prepareCurveData( edges, histdata );
% for xx = 1:20 
%     if yData(xx,1) >= abs(100)
%         yData(xx,1) = 0;
%     end
% end
% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [24.4274436090225 3790000 857734.333125455];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
G = figure( 'Name', 'Gauss1FitDt','visible','off' );
yFitted = feval(fitresult,xData);
[pks,loc] = findpeaks(yFitted,xData);
plot( fitresult,xData,yData);
hold on
plot(loc,pks,'vr','MarkerEdgeColor','r','MarkerFaceColor','r');
axis([0 max(edges) 0 (max(yFitted)+10)])
legend( Name,'Gaussian Fit',['Max: ' num2str(loc) ' nm²/s'] ,'Location', 'NorthEast' );
legend ('boxoff')
% Label axes
xlabel 'D_t [nm²/s]'
ylabel 'Particlecounts'
set(gcf,'color','w');
grid on
locations = loc;
export_fig([output '_Gauss2Fit'],G,'-pdf','-tiff','-r300','-painters')
ExpData = array2table([xData,yData,yFitted]);
headers = {'Diameter' 'MeasuredData' 'FittedData'};
ExpData.Properties.VariableNames = headers;
writetable(ExpData,[output '_Gauss2Fit']);
close (G);

end

