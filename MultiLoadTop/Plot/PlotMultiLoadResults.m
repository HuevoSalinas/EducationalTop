function [] = PlotMultiLoadResults(method, MaxIter, hist_J, hist_V, hist_w)

switch method
    case 'pareto'
        Nl = size(hist_J,3);
        H = linspace(0, 1, Nl + 1);
        H = H(1: Nl)';
        palette = hsv2rgb([H ones(Nl, 1) ones(Nl,1)]);
        [~, ind] = min(max(hist_J(:,MaxIter,:), [], 3));
        hist_J_plot = reshape(hist_J(ind,:,:), [MaxIter+1, Nl]);
        hist_V_plot = reshape(hist_V(ind,:),[MaxIter+1, 1]) ;
        hist_w_plot = [];

        % Pareto front
        if size(hist_w,2) == 2
            J1 = hist_J(:, MaxIter, 1);
            J2 = hist_J(:, MaxIter, 2);
            alpha = hist_w(:, 1);
            minJ = min(max(hist_J(:,MaxIter,:),[],3));
            figure, hold on, grid on, box on
            plot(alpha, J1,'LineWidth',1,'Color',[0.64 0.29 0.64])
            plot(alpha, J2,'LineWidth',1,'Color',[1 0.29 0])
            plot(alpha, alpha.*J1 + (1-alpha).*J2,'LineWidth',1,'Color',[0 0 0])
            
            legend('$J_1(w)$', '$J_2(w)$', '$ \theta(w)$')
            ylim([0, minJ*1.5])
            set(gca,'YScale','log'), xlabel('$w_1$'), ylabel('objective J')
            ApplyLatexInterpreter(gcf)
        
        elseif size(hist_w,2) == 3
            PlotParetoFront3D(hist_J(:, MaxIter, :), hist_w, 'map',true)
        end
            
    case 'pareto-search'
        Nl = size(hist_J,3);
        H = linspace(0, 1, Nl + 1);
        H = H(1: Nl)';
        palette = hsv2rgb([H ones(Nl, 1) ones(Nl,1)]);
        [~, ind] = min(max(hist_J(:,MaxIter,:), [], 3));
        hist_J_plot = reshape(hist_J(ind,:,:), [MaxIter+1, Nl]);
        hist_V_plot = reshape(hist_V(ind,:),[MaxIter+1, 1]) ;
        hist_w_plot = [];
 
        % Pareto front
        if size(hist_w,2) == 3
            PlotParetoFront3D(hist_J(:, MaxIter, :), hist_w, 'trajectory', true)
        end
    case 'bounded-dual'
        Nl = size(hist_J,2);
        H = linspace(0, 1, Nl + 1);
        H = H(1: Nl)';
        palette = hsv2rgb([H ones(Nl, 1) ones(Nl,1)]);
        hist_J_plot = hist_J;
        hist_V_plot = hist_V;
        hist_w_plot = hist_w;
       
        % Pareto front
        if size(hist_w,2) == 3
            PlotParetoFront3D(hist_J(1:MaxIter,:), hist_w(1:MaxIter,:), 'trajectory', false)
        end

    otherwise
        Nl = size(hist_J,2);
        H = linspace(0, 1, Nl + 1);
        H = H(1: Nl)';
        palette = hsv2rgb([H ones(Nl, 1) ones(Nl,1)]);
        hist_J_plot = hist_J;
        hist_V_plot = hist_V;
        hist_w_plot = hist_w;     
end

 % Objective function 
figure, hold on, grid on, box on
for i = 1:size(hist_J_plot, 2)
    plot(0:MaxIter,hist_J_plot(:, i),'LineWidth',1, ...
        'Color', palette(i,:), 'DisplayName', sprintf('Obj %i', i) )
end
legend
set(gca,'YScale','log'), xlabel('iteration'), ylabel('objective J')
ApplyLatexInterpreter(gcf)

figure, hold on, grid on, box on
plot(0:MaxIter,max(hist_J_plot, [], 2),'LineWidth',1)
set(gca,'YScale','log'), xlabel('iteration'), ylabel('max objective J')
ApplyLatexInterpreter(gcf)

% Weights
if ~isempty(hist_w_plot)
    figure, hold on, grid on, box on
    for i = 1:size(hist_w_plot, 2)
        plot(0:MaxIter,hist_w_plot(:, i),'LineWidth',1, ...
            'Color', palette(i,:), 'DisplayName', sprintf('Obj %i', i) )
    end
    legend
    xlabel('iteration'), ylabel('weight')
    ApplyLatexInterpreter(gcf)
end

% Volume fraction
figure, hold on, grid on, box on
plot(0:MaxIter,hist_V_plot,'LineWidth',1,'Color',[1.00 0.50 0.10])
xlabel('iteration'), ylabel('volume fraction')
ApplyLatexInterpreter(gcf)
end

function [] = PlotParetoFront3D(hist_J, hist_w, plot_type, color_on)
TRIwidth = 1/sin(pi/3);
ticks = 0:0.2:1;
ticksize = 1/20;
CBscale = 0.85;
GCAxoffset = -0.15;
J1 = hist_J(:, 1);
J2 = hist_J(:, 2);
J3 = hist_J(:, 3);
w1 = hist_w(:,1);
w2 = hist_w(:,2);
w3 = hist_w(:,3);
TRI = [-2/3*cos(pi/6) -2/3*sin(pi/6);
        2/3*cos(pi/6) -2/3*sin(pi/6);
        0             2/3      ];

H = max([J1, J2, J3], [], 2);
theta = w1.*J1 + w2.*J2 + w3.*J3;
capH = min(H) * 1.5;
captheta = max(theta) * 0.6;
H = min(H, capH);
theta = max(theta, captheta);
xi = [w2(:) w3(:)];
NODE = nan(size(xi,1),2);
for i=1:size(xi,1)
    NODE(i,:) = ShapeFuncT3(xi(i,:))*TRI;
end
switch plot_type
    case 'map'
        ELEM = delaunay(NODE);
    case 'trajectory'
        xi_map = [1 0;
              0 1;
              0 0];
        NODE_map = nan(size(xi_map,1),2);
        for i=1:size(xi_map,1)
            NODE_map(i,:) = ShapeFuncT3(xi_map(i,:))*TRI;
        end
        ELEM_map = [1 2 3];
    otherwise
         error('Unknown plot type.')
end

figure, hold on, axis equal, grid on, box on
switch plot_type
    case 'map'
        patch('Vertices',NODE,'Faces',ELEM,'EdgeColor','none',...
            'FaceColor','interp','FaceVertexCData',H)
        [~,ind]= min(H);
    case 'trajectory'
        patch('Vertices',NODE_map,'Faces',ELEM_map,'EdgeColor','none',...
                'FaceColor',[0.9 0.9 0.9])
        p = plot(NODE(:,1),NODE(:,2), LineWidth=2);
        if color_on
            nc = size(xi,1);
            color_cuts = linspace(min(H), max(H), 16);
            color_inds = ones(nc,1);
            for i = 1:16
                color_inds(H>=color_cuts(i)) = i;
            end
            cm = sky(16);
            colors = [uint8(cm(color_inds,:)*255) uint8(ones(nc,1))].';
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colors)
        end
        ind = size(xi,1);

    otherwise
         error('Unknown plot type.')
end

DrawWGrid(TRIwidth,ticks,ticksize)
plot(NODE(ind,1),NODE(ind,2),'kx')
if NODE(ind,1)<=0
    aux = sprintf('$\\;\\;\\leftarrow$ optimum');
    text(NODE(ind,1),NODE(ind,2),aux)
else
    aux = sprintf('optimum $\\rightarrow\\;\\;$');
    text(NODE(ind,1),NODE(ind,2),aux,'HorizontalAlignment','right')
end
axis off, axis tight, axis(1.25*[-TRIwidth/2 TRIwidth/2 -1/3 2/3])
ApplyLatexInterpreter(gcf)

if color_on
    colormap(sky(16))
    CBhand = colorbar;
    aux = CBhand.Position;
    aux(2) = aux(2) + 0.333*(1-CBscale)*aux(4);
    aux(4) = CBscale*aux(4);
    CBhand.Position = aux;
    aux = gca;
    aux.Position(1) = aux.Position(1) + GCAxoffset;
    old_label = CBhand.TickLabels{end};
    CBhand.TickLabels{end} = sprintf('above %s',old_label);
    CBhand.Title.String='$ max\{J_i(w)\} $';
    set(CBhand.XLabel,{'String','Rotation','Position'},{'XLabel',0,[0.5 -0.01]});
end
ApplyLatexInterpreter(gcf)


figure, hold on, axis equal, grid on, box on
switch plot_type
    case 'map'
        patch('Vertices',NODE,'Faces',ELEM,'EdgeColor','none',...
            'FaceColor','interp','FaceVertexCData',theta)
        [~,ind]= max(theta);
    case 'trajectory'
        patch('Vertices',NODE_map,'Faces',ELEM_map,'EdgeColor','none',...
                'FaceColor',[0.9 0.9 0.9])
        p = plot(NODE(:,1),NODE(:,2), LineWidth=2);
        if color_on
            nc = size(xi,1);
            color_cuts = linspace(min(theta), max(theta), 16);
            color_inds = ones(nc,1);
            for i = 1:16
                color_inds(theta>=color_cuts(i)) = i;
            end
            cm = summer(16);
            colors = [uint8(cm(color_inds,:)*255) uint8(ones(nc,1))].';
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',colors)
         end
         ind = size(xi,1);
    otherwise
         error('Unknown plot type.')
end

DrawWGrid(TRIwidth,ticks,ticksize)
plot(NODE(ind,1),NODE(ind,2),'kx')
if NODE(ind,1)<=0
    aux = sprintf('$\\;\\;\\leftarrow$ optimum');
    text(NODE(ind,1),NODE(ind,2),aux)
else
    aux = sprintf('optimum $\\rightarrow\\;\\;$');
    text(NODE(ind,1),NODE(ind,2),aux,'HorizontalAlignment','right')
end
axis off, axis tight, axis(1.25*[-TRIwidth/2 TRIwidth/2 -1/3 2/3])
ApplyLatexInterpreter(gcf)

if color_on
    colormap(summer(16))
    CBhand = colorbar;
    aux = CBhand.Position;
    aux(2) = aux(2) + 0.333*(1-CBscale)*aux(4);
    aux(4) = CBscale*aux(4);
    CBhand.Position = aux;
    aux = gca;
    aux.Position(1) = aux.Position(1) + GCAxoffset;
    old_label = CBhand.TickLabels{1};
    CBhand.TickLabels{1} = sprintf('below %s',old_label);
    CBhand.Title.String='$ \theta(w) $';
    set(CBhand.XLabel,{'String','Rotation','Position'},{'XLabel',0,[0.5 -0.01]});
end
ApplyLatexInterpreter(gcf)



end

function [] = DrawWGrid(TRIwidth,ticks,ticksize)

% W lines
for i=1:length(ticks)
    my_width = (1-ticks(i))*TRIwidth;
    U = my_width/2*[-1 1] + ticksize*[0 1];
    V = (ticks(i)-1/3)*[1 1];
    plot(U,V,'k')
    aux = sprintf('$\\;\\;w_3=%g$',ticks(i));
    text(U(end),V(end),aux,'Interpreter','latex')

    U = my_width/2*[-1 1] + ticksize*[0 1];
    V = (ticks(i)-1/3)*[1 1];
    X = [cos(-pi*4/3) sin(-pi*4/3)] * vertcat(U,V);
    Y = [-sin(-pi*4/3) cos(-pi*4/3)] * vertcat(U,V);
    plot(X,Y,'k')
    aux = sprintf('$w_2=%g\\;\\;$',ticks(i));
    text(X(end),Y(end),aux,...
        'HorizontalAlignment','right','Rotation',60)

    U = my_width/2*[-1 1] + ticksize*[0 1];
    V = (ticks(i)-1/3)*[1 1];
    X = [cos(pi*4/3) sin(pi*4/3)] * vertcat(U,V);
    Y = [-sin(pi*4/3) cos(pi*4/3)] * vertcat(U,V);
    plot(X,Y,'k')
    aux = sprintf('$w_1=%g\\;\\;$',ticks(i));
    text(X(end),Y(end),aux,...
        'HorizontalAlignment','right','Rotation',-60)
end

end
