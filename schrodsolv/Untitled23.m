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
    T_all{i} = data{1}(end);
    R_all{i} = data{2}(end);
end;
T_array = cell2mat(T_all);
R_array = cell2mat(R_all);
S = T_array + R_array;

plot(S)
