function [optoutput] = PlotCurrentDensity(NODE, ELEM, x, optoutput, iter)

FaceColor = ones(length(x), 1);
FaceAlpha = x;

if (iter==1)
    optoutput.topfig = figure('Color','w');
    axis equal, hold on, grid on, box on, axis off
    optoutput.toppatch = PlotDensity(NODE,ELEM,FaceColor,FaceAlpha);
    PlotWindowLimits(NODE)
    colormap(flipud(gray)), set(gca,'CLim',[0 1])
    
    optoutput.topfigtitle = title(sprintf('iteration = %03.0f',iter-1));
    ApplyLatexInterpreter(gcf)
else
    set(optoutput.toppatch,'FaceVertexAlphaData',x)
    set(optoutput.topfigtitle,'String',sprintf('iteration = %03.0f',iter-1))
end
   
        
return 