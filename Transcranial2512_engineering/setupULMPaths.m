function [dataDir,localPath,trackspath,savingpath,fileList,allFileNames,nBuffers] = setupULMPaths(baseDir, expName, runName)
%SETULMPATHS 统一设置项目路径并创建输出目录。
%
%   [dataDir,localPath,trackspath,savingpath,fileList,allFileNames,nBuffers] = setupULMPaths(...)
%   只封装路径相关逻辑。
%   采集参数和算法参数仍保留在 main 脚本中，便于直接检查和调参。

if nargin < 1 || isempty(baseDir)
    baseDir = fileparts(mfilename('fullpath'));
end

if nargin < 2 || isempty(expName)
    expName = 'default';
end

if nargin < 3 || isempty(runName)
    runName = expName;
end

workDir = fullfile(baseDir, 'SR');
dataDir = fullfile(baseDir, 'data', expName);
localPath = fullfile(workDir, 'Locals', runName);
trackspath = fullfile(workDir, 'Tracks', runName);
savingpath = fullfile(workDir, 'Mats', runName);

if ~exist(dataDir, 'dir')
    originalDataDir = fullfile(fileparts(baseDir), 'Transcranial2512', 'data', expName);
    if exist(originalDataDir, 'dir')
        % 工程副本中不重复保存大数据；若本地 data 不存在，则读取原始项目的数据目录。
        dataDir = originalDataDir;
    end
end

% 将项目根目录及子目录加入 MATLAB 搜索路径，保证算法函数可被调用。
addpath(genpath(baseDir));

if ~exist(dataDir, 'dir')
    error('Input data folder not found: %s\nPlease check the data path.', dataDir);
end

filePattern = fullfile(dataDir, '*.mat');
fileList = dir(filePattern);
if isempty(fileList)
    error('No .mat files found in data folder: %s', dataDir);
end
allFileNames = {fileList.name};
nBuffers = length(allFileNames);

outputPaths = {localPath, trackspath, savingpath};
for iPath = 1:numel(outputPaths)
    if ~exist(outputPaths{iPath}, 'dir')
        % 输出目录不存在时自动创建，避免每个 main 脚本重复写 mkdir 逻辑。
        mkdir(outputPaths{iPath});
        fprintf('>>> Created folder: %s\n', outputPaths{iPath});
    end
end

end
