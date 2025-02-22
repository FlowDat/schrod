clear all
clc
% Get list of files that match the pattern
filePattern = 'RT-0.1-*.txt';
files = dir(filePattern);

% Initialize arrays to store x and y values
x_values = [];
y_values = [];

% Loop through each file
for i = 1:length(files)
    % Extract the number at the end of the filename
    fileName = files(i).name;
    number = str2double(regexp(fileName, '\d+', 'match', 'once'));
    x_values = [x_values, number];
    
    % Read the file content
    filePath = fullfile(folder, fileName);
    fileData = importdata(filePath);
    
    % Sum the last two lines of data
    lastTwoSum = sum(fileData(end-1:end));
    y_values = [y_values, lastTwoSum];
end

% Plot the data
plot(x_values, y_values, 'o-', 'LineWidth', 1.5);
xlabel('X Value');
ylabel('Sum of Last Two Lines');
title('Sum of Data in Last Two Lines vs X Value');
grid on;
