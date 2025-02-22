% Directory containing the files
filePattern = 'RT-*.txt';
files = dir(filePattern);

% Initialize cell arrays to store data
numFiles = length(files);
T_all = cell(numFiles, 1);
R_all = cell(numFiles, 1);

% Define different plot styles
plotStyles = {'-o', '-s', '-d', '--o', '--s', '--d', ':o', ':s', ':d', '-.*'};
numStyles = length(plotStyles);

% Loop through each file
for i = 1:numFiles
    % Get the file name
    fileName = files(i).name;
    
    % Open the file for reading
    fileID = fopen(fileName, 'r');
    
    % Read the header line (optional)
    header = fgetl(fileID);
    
    % Read the data from the file
    data = textscan(fileID, '%f %f', 'Delimiter', '\t');
    
    % Close the file
    fclose(fileID);
    
    % Store T and R in cell arrays
    T_all{i} = data{1};
    R_all{i} = data{2};
end
 
% Plot the data
figure;
hold on;
tt=1:19;
for i = 1:numFiles
    % Cycle through plot styles
    plotStyle = plotStyles{mod(i-1, numStyles) + 1};
    plot(tt,T_all{i} + R_all{i}, plotStyle,'LineWidth', 1.5, 'DisplayName', sprintf('File %s', files(i).name));
    %plot(tt,T_all{i}, plotStyle,'LineWidth', 1.5, 'DisplayName', sprintf('File %s', files(i).name));
end
%ylim([0.166, 0.1698]);
hold off;
xlabel('time');
ylabel('Probability');
%title('Plot of Survival Probabilities');
legend show;
grid on;
