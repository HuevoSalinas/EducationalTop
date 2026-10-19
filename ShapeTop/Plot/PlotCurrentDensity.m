function [optoutput] = PlotCurrentDensity(NODE, ELEM, x, optoutput, iter)

Ne = size(ELEM,1);
FaceColor = repmat([0 0 0],Ne,1);
FaceAlpha = x;

if (iter==1)
    optoutput.topfig = figure('Color','w');
    axis equal, hold on, grid on, box on, axis off
    optoutput.toppatch = PlotDensity(NODE,ELEM,FaceColor,FaceAlpha);
    edges = [ELEM(:,[1 2]); ELEM(:,[2 3]); ELEM(:,[3 4]); ELEM(:,[4 1])];
    [uniqueEdges,~,edgeID] = unique(sort(edges,2),'rows');
    edgeCount = accumarray(edgeID,1);
    boundary = uniqueEdges(edgeCount==1,:);
    plot([NODE(boundary(:,1),1) NODE(boundary(:,2),1)]',...
         [NODE(boundary(:,1),2) NODE(boundary(:,2),2)]',...
         'Color',[0.11 0.46 0.74],'LineWidth',0.5)
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