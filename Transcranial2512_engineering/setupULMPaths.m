function [dataDir,localPath,trackspath,savingpath,fileList,allFileNames,nBuffers] = setupULMPaths(baseDir, expName)
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

workDir = fullfile(baseDir, 'SR');
dataDir = fullfile(baseDir, 'data', expName);
localPath = fullfile(workDir, 'Locals', expName);
trackspath = fullfile(workDir, 'Tracks', expName);
savingpath = fullfile(workDir, 'Mats', expName);

% 将项目根目录及子目录加入 MATLAB 搜索路径，保证算法函数可被调用。
addpath(genpath(baseDir));

if ~exist(dataDir, 'dir')
    error('【错误】找不到输入文件夹：%s\n请检查数据位置！', dataDir);
end

filePattern = fullfile(dataDir, '*.mat');
fileList = dir(filePattern);
if isempty(fileList)
    error('【错误】在路径 [%s] 下未找到 .mat 文件！', dataDir);
end
allFileNames = {fileList.name};
nBuffers = length(allFileNames);

outputPaths = {localPath, trackspath, savingpath};
for iPath = 1:numel(outputPaths)
    if ~exist(outputPaths{iPath}, 'dir')
        % 输出目录不存在时自动创建，避免每个 main 脚本重复写 mkdir 逻辑。
        mkdir(outputPaths{iPath});
        fprintf('>>> 已创建文件夹: %s\n', outputPaths{iPath});
    end
end

end
