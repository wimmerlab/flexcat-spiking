function figsave (handle, filename, figformat)
%FIGSAVE  Save figure in the specified formats
%   FIGSAVE(handle, filename, figformat) save the figure in fig and/or png
%   format or don't save anything if 'figformat' is empty

for i=1:numel(figformat)
    switch (figformat{i})
        case 'png'
            saveas(handle, [filename '.png'], 'png');
        case 'fig'
            saveas(handle, [filename '.fig'], 'fig');
        otherwise
            error('unknown format')
    end
end
