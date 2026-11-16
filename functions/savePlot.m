function savePlot(obj, fileName)
    % h = gcf;
    % set(h, 'PaperType', 'uslegal');
    % set(h, 'PaperPositionMode', 'auto');
    % set(h, 'PaperOrientation', 'landscape');
    % set(h, 'Position', [10 60 1600 900]);
    % fullFileName = fullfile(obj.folder, '\results\', fileName);
    % print(fullfile(fullFileName), '-dpng', '-r300');

    % decide where to put it
    outDir = fullfile(obj.folder, 'results');
    % create it if needed
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % build full filename and save
    outFile = fullfile(outDir, [fileName, '.png']);
    print(outFile, '-dpng', '-r300');    
end