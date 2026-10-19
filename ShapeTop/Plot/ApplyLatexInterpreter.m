function ApplyLatexInterpreter(fig)
    
% Get all axes in the figure
axs = findall(fig, 'Type', 'axes');

for ax = axs'
    % Apply LaTeX interpreter to axis labels and title
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.ZLabel.Interpreter = 'latex';
    ax.Title.Interpreter = 'latex';

    % Modify tick labels using LaTeX
    ax.TickLabelInterpreter = 'latex';

    % Apply LaTeX interpreter to legend
    if ~isempty(ax.Legend)
       ax.Legend.Interpreter = 'latex';
    end

    % Apply LaTeX interpreter to colorbars
    cbar = findall(fig, 'Type', 'colorbar', 'Axes', ax);
    for cb = cbar(:)'
        cb.Title.Interpreter = 'latex';
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';
    end

    % Apply LaTeX interpreter to texts
    textElements = findall(ax, 'Type', 'text'); % Find all text objects
    for txt = textElements'
        txt.Interpreter = 'latex';
    end
end

end
