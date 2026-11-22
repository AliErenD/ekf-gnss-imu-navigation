function save_fig(filename)
%SAVE_FIG Save current figure as PNG reliably in -batch (headless) mode.
set(gcf, 'Renderer', 'painters');
print(gcf, filename, '-dpng', '-r150');
end
