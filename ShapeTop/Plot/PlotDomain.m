function [handle] = PlotDomain(NODE,ELEM,SUPP,LOAD,FaceColor,FaceAlpha)
narginchk(2,6)

figure, axis equal, hold on, grid on, box on
if ((nargin<5) || isempty(FaceColor))
    FaceColor = [0.9 0.9 1];
end
if (nargin<6)
    FaceAlpha = 1;
end
handle = patch('Faces',ELEM,'Vertices',NODE,'FaceColor','flat',...
               'FaceVertexCData',FaceColor,'FaceAlpha',FaceAlpha);

if (nargin>=3 && ~isempty(SUPP))
    plot(NODE(SUPP(:,1),1),NODE(SUPP(:,1),2),'k>','MarkerSize',8);
end
if (nargin>=4 && ~isempty(LOAD))
    LOADval = vecnorm(LOAD(:,2:end),2,2);
    LOADval = LOADval / max(LOADval);
    for i=1:length(LOADval)
        if LOADval(i)>0
            plot(NODE(LOAD(i,1),1),NODE(LOAD(i,1),2),'r^','MarkerSize',15*LOADval(i));
        end
    end
end
PlotWindowLimits(NODE)
ApplyLatexInterpreter(gcf)
drawnow

end
