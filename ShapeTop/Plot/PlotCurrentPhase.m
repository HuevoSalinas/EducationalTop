function [optoutput] = PlotCurrentPhase(NODE, ELEM, x, opt, optoutput, iter)

PhaseColor = repmat(0.5,length(x),1);
PhaseColor(opt.flagPureVoid,:) = 0;
PhaseColor(opt.flagPureSolid,:) = 1;
PhaseAlpha = ones(length(x),1);
PhaseAlpha(opt.flagPureVoid,:) = 0;

if (iter==1)
   optoutput.phasefig = figure('Color','w');
   axis equal, hold on, grid on, box on, axis off
   PlotWindowLimits(NODE)
   optoutput.phasepatch = PlotDensity(NODE,ELEM,PhaseColor,PhaseAlpha);
    edges = [ELEM(:,[1 2]); ELEM(:,[2 3]); ELEM(:,[3 4]); ELEM(:,[4 1])];
    [uniqueEdges,~,edgeID] = unique(sort(edges,2),'rows');
    edgeCount = accumarray(edgeID,1);
    boundary = uniqueEdges(edgeCount==1,:);
    plot([NODE(boundary(:,1),1) NODE(boundary(:,2),1)]',...
         [NODE(boundary(:,1),2) NODE(boundary(:,2),2)]',...
         'Color',[0.11 0.46 0.74],'LineWidth',0.5)
   colormap([1 1 1; 0.5 0.5 0.5; 0 0 0])

   ColorbarHand = colorbar;
   set(ColorbarHand,'Ticks',[1/6 1/2 5/6],'TickLabels',{'void','trans','solid'})
   fig1posit = get(optoutput.phasefig,'Position');
   set(optoutput.phasefig,'Position',[fig1posit(1)-fig1posit(3) fig1posit(2:4)])
   aux = sprintf('Material Phase: $\\mu = %.2f$',opt.current_mu);
   optoutput.phasefigtitle = title(aux);
   ApplyLatexInterpreter(gcf)
   drawnow
else
    set(optoutput.phasepatch,'FaceVertexCData',PhaseColor,'FaceVertexAlphaData',PhaseAlpha)
    aux = sprintf('Material Phase: $\\mu = %.2f$',opt.current_mu);
    set(optoutput.phasefigtitle,'String',aux);
    drawnow
end
        
end