function [optoutput] = PlotCurrentDensity(NODE, ELEM, x, optoutput, iter)

Ne = size(ELEM,1);
FaceColor = repmat([0 0 0],Ne,1);
FaceAlpha = x;

if (iter==1)
    optoutput.topfig = figure('Color','w');
    axis equal, hold on, grid on, box on, axis off
    optoutput.toppatch = PlotDensity(NODE,ELEM,FaceColor,FaceAlpha);
    PlotWindowLimits(NODE)
    
    optoutput.topfigtitle = title(sprintf('iteration = %03.0f',iter-1));
    ApplyLatexInterpreter(gcf)
    drawnow
else
    set(optoutput.toppatch,'FaceVertexAlphaData',FaceAlpha)
    set(optoutput.topfigtitle,'String',sprintf('iteration = %03.0f',iter-1))
    drawnow
end
   
        
end