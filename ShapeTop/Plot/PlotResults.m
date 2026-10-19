function [] = PlotResults(MaxIter, hist_J, hist_V)

% Objective function 
figure, hold on, grid on, box on
plot(0:MaxIter,hist_J,'LineWidth',1,'Color',[0.64 0.29 0.64])
set(gca,'YScale','log'), xlabel('iteration'), ylabel('objective J')
ApplyLatexInterpreter(gcf)

% Volume fraction
figure, hold on, grid on, box on
plot(0:MaxIter,hist_V,'LineWidth',1,'Color',[1.00 0.50 0.10])
xlabel('iteration'), ylabel('volume fraction')
ApplyLatexInterpreter(gcf)

end