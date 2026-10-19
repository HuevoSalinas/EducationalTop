function [] = PlotDisplacementsDiffTOP(NODE,ELEM,D1,x1,D2,x2,dispFactor)

values1 = vecnorm(D1,2,2);
values2 = vecnorm(D2,2,2);
values = [values1 values2];
x = {x1,x2};
DNODE1 = NODE + dispFactor*D1;
DNODE2 = NODE + dispFactor*D2;
DNODE = {DNODE1,DNODE2};
Ne = size(ELEM,1);
NpE = size(ELEM,2);
NODEx1 = accumarray(ELEM(:),repmat(x1,NpE,1)) ./ ...
          accumarray(ELEM(:),ones(Ne*NpE,1));
NODEx2 = accumarray(ELEM(:),repmat(x2,NpE,1)) ./ ...
          accumarray(ELEM(:),ones(Ne*NpE,1));
NODEx = {NODEx1,NODEx2};
maxval = max([values1(NODEx1>0.5); values2(NODEx2>0.5)]);
if isempty(maxval) || ~isfinite(maxval) || maxval<=0
    maxval = 1;
end
limit_nodes = [NODE; DNODE1(NODEx1>0.5,:); DNODE2(NODEx2>0.5,:)];

figure
for i=1:2
    subplot(2,1,i)
    hold on, axis equal, grid off, box on
    patch('Faces',ELEM,'Vertices',DNODE{i},...
          'FaceColor','interp','FaceVertexCData',values(:,i),...
          'EdgeColor','none','FaceAlpha','flat','FaceVertexAlphaData',x{i});
    title(sprintf('Maximum displacement = %g',max(values(:,i))))
    clim([0 maxval])
    colormap(plasma)
    colorbar
    PlotWindowLimits(limit_nodes,0)
end
ApplyLatexInterpreter(gcf)
drawnow

end
