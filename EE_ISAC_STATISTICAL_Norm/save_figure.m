function saveAllFigures()
    % saveAllFigures Finds all open figures and saves them in .png and .fig
    % formats inside a timestamped subfolder within results/figures/.

    % 1. Create the timestamped folder
    % Format: YYYY-MM-DD_HH-MM-SS (e.g., 2026-04-13_11-00-11)
    timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
  % 1. Get the parent directory and the current folder name
% If pwd is 'D:\EE_ISAC_BoundedCSI_Norm':
%   parentDir becomes 'D:\'
%   currentFolder becomes 'EE_ISAC_BoundedCSI_Norm'
[parentDir, currentFolder, ~] = fileparts(pwd);

% 2. Build the base directory path using the EXACT source folder name
% Navigates to D:\Results\EE_ISAC_BoundedCSI_Norm
baseDir = fullfile(parentDir, 'Results', currentFolder);

% 3. Add the timestamp for this specific run
saveDir = fullfile(baseDir, timestamp);
    % Create the folder if it doesn't exist
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    % 2. Find all currently open figures
    figs = findall(0, 'Type', 'figure');
    
    if isempty(figs)
        disp('No open figures found to save.');
        return;
    end
    
    fprintf('Saving %d figures to: %s\n', length(figs), saveDir);

    % 3. Loop through each figure and save it
    for i = 1:length(figs)
        f = figs(i);
        
        % Generate a smart filename
        if isempty(f.Name)
            % If the figure has no title, use its figure number
            fileName = sprintf('Figure_%d', f.Number);
        else
            % If it has a name, sanitize it to be filesystem-safe
            fileName = regexprep(f.Name, '[^\w-]', '_');
        end
        
        % Define the full paths for both formats
        pngPath = fullfile(saveDir, [fileName, '.png']);
        figPath = fullfile(saveDir, [fileName, '.fig']);
        
        % Save as high-resolution PNG
        exportgraphics(f, pngPath, 'Resolution', 300);
        
        % Save as editable MATLAB Figure
        savefig(f, figPath);
    end
    
    disp('All figures saved successfully.');
end
