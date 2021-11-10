function [] = OpenMostRecentDate()
% This function reads most recent date.

%% 
condintionName = 'ConeIsolating';
testFileName = sprintf('testImageDataCheck_%s',conditionName);
testFileDir = dir(append(testFileName,'*'));
allFiles = zeros(length(testFileDir));
 
for ii = 1:length(testFileDir)
    fileDirDate_temp = testFileDir(ii).date;
    allFiles(ii) =datenum(fileDirDate_temp);
end

[tmp i]=max(allFiles);
    
load(allFiles(i).name)

end