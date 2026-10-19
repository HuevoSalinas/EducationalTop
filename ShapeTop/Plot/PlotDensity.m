function [handle] = PlotDensity(NODE,ELEM,FaceColor,FaceAlpha)

handle = patch('Faces',ELEM,'Vertices',NODE,'FaceColor','flat', ...
               'FaceVertexCData',FaceColor, 'FaceAlpha','flat', ...
               'AlphaDataMapping','none','FaceVertexAlphaData', FaceAlpha,...
               'EdgeColor','none');

drawnow

end
