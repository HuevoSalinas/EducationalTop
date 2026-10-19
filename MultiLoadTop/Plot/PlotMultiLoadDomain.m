function [handle] = PlotMultiLoadDomain(NODE,ELEM,SUPP,LOAD,FaceColor,FaceAlpha)
narginchk(2,6)

figure, axis equal, hold on, grid on, box on
if ((nargin<5) || isempty(FaceColor))
    FaceColor = [0.9 0.9 1];
end
if (nargin<6)
    FaceAlpha = 1;
end

handle = patch('Faces',ELEM,'Vertices',NODE,'FaceColor','flat','FaceVertexCData',FaceColor,'FaceAlpha',FaceAlpha);

if (nargin>=3 && ~isempty(SUPP))
    if size(NODE,2)==2
        plot(NODE(SUPP(:,1),1),NODE(SUPP(:,1),2),'k>','MarkerSize',8);   
    else
        error('Function can only take 2D nodes.')
    end
end

if (nargin>=4 && ~isempty(LOAD))
    Nl = size(LOAD,2);
    H = linspace(0, 1, Nl + 1);
    H = H(1: Nl)';
    palette = hsv2rgb([H ones(Nl, 1) ones(Nl,1)]);
    for e = 1:Nl
        LOADval = vecnorm(LOAD{e}(:,2:end),2,2);
        LOADval = LOADval / max(LOADval);
        if (~isempty(LOAD{e}))
            if size(NODE,2)==2
                for i=1:length(LOADval)
                    if (LOADval(i)>0)
                        plot(NODE(LOAD{e}(i,1),1),NODE(LOAD{e}(i,1),2),'^', ...
                            'Color', palette(e, :), 'MarkerSize',15*LOADval(i));
                    end
                end
            else
                error('Function can only take 2D nodes.')
            end
        end
        
    end
    colormap(palette)
    ColorbarHand = colorbar;
    set(ColorbarHand,'Ticks',(1/(2*Nl)):1/(Nl):1,'TickLabels',["Load " + string(1:Nl)])
end

PlotWindowLimits(NODE)
ApplyLatexInterpreter(gcf)
drawnow
end
