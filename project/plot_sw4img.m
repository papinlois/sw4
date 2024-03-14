disp('Need to figure what each files/plots are meant to represent')
disp('so we can choose what we want to see')

% Specify the directory containing the image files
directory = '/Users/lpapin/Documents/sw4-master/mine/vci-results/';

% Get a list of all files in the directory
files = dir(directory);

% Loop through each file in the directory
for i = 4:length(files)
    % Get the file name
    filename = files(i).name;
    
    % Construct the full file path
    filepath = fullfile(directory, filename);

    % Read the image using the readimage function
    [im, x, y, z, plane, t, timestring, npatches] = readimage(filepath, 1, 0);

    % Plot the image
    figure;
    if plane == 0
        imagesc(y, z, im);
        colormap(parula); % Set colormap
        colorbar; % Display colorbar
        xlabel('Y');
        ylabel('Z');
        title(['SW4 Image: ' filename]);
    elseif plane == 1
        imagesc(x, z, im);
        colormap(parula); % Set colormap
        colorbar; % Display colorbar
        xlabel('X');
        ylabel('Z');
        title(['SW4 Image: ' filename]);
    else
        imagesc(x, y, im);
        colormap(parula); % Set colormap
        colorbar; % Display colorbar
        xlabel('X');
        ylabel('Y');
        title(['SW4 Image: ' filename]);
    end
end
