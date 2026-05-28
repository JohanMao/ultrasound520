function viewULMLocalization(blockIdx, varargin)
%VIEWULMLOCALIZATION 显示指定 block 内定位点的空间分布。
%
%   viewULMLocalization(blockIdx)
%   显示指定 block 内所有定位点。
%
%   viewULMLocalization(blockIdx, 'FrameRange', [1 20])
%   只显示指定 block 内第 1 到 20 帧的定位点。
%
%   viewULMLocalization(blockIdx, 'FrameStart', 101, 'FrameWindow', 50)
%   只显示指定 block 内第 101 到 150 帧的定位点。
%
%   viewULMLocalization(blockIdx, 'FrameWindow', 20)
%   从该 block 的第一帧开始，只显示 20 帧定位点。
%
%   viewULMLocalization(..., 'ExpName', '20250411flow')
%   指定实验数据文件夹名称。
%
%   viewULMLocalization(..., 'RunName', '20250411flow_svd10_np200')
%   指定要读取的 SR/Locals 输出结果文件夹名称。

if nargin < 1 || isempty(blockIdx)
    blockIdx = 1;
end

parser = inputParser;
parser.addParameter('ExpName', '20250411flow', @(x) ischar(x) || isstring(x));
parser.addParameter('RunName', '', @(x) ischar(x) || isstring(x));
parser.addParameter('PointSize', 18, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('FrameRange', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
parser.addParameter('FrameStart', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
parser.addParameter('FrameWindow', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
parser.addParameter('Mode', 'scatter', @(x) any(strcmpi(char(x), {'scatter', 'density'})));
parser.parse(varargin{:});

expName = char(parser.Results.ExpName);
runName = char(parser.Results.RunName);
pointSize = parser.Results.PointSize;
frameRange = parser.Results.FrameRange;
frameStart = parser.Results.FrameStart;
frameWindow = parser.Results.FrameWindow;
displayMode = lower(char(parser.Results.Mode));

if isempty(runName)
    runName = expName;
end

baseDir = fileparts(mfilename('fullpath'));
[dataDir, localPath, ~, ~, ~, allFileNames, nBuffers] = setupULMPaths(baseDir, expName, runName);

if blockIdx < 1 || blockIdx > nBuffers
    error('blockIdx is out of range. Available block count: %d.', nBuffers);
end

localFile = fullfile(localPath, sprintf('Locals_%03d.mat', blockIdx));

if ~exist(localFile, 'file')
    error('Localization result file not found: %s\nPlease run ULM_main2D.m first to generate Locals files.', localFile);
end

localStruct = load(localFile, 'MatTracking');
if ~isfield(localStruct, 'MatTracking')
    error('Localization result file does not contain variable MatTracking: %s', localFile);
end

MatTracking = localStruct.MatTracking;
if isempty(MatTracking)
    warning('Localization result is empty for block %03d.', blockIdx);
end

if size(MatTracking, 2) < 4
    error('MatTracking must contain at least four columns: intensity, axial, lateral, frame.');
end

availableFrames = unique(MatTracking(:, 4));
if isempty(availableFrames)
    minFrame = 1;
    maxFrame = 1;
else
    minFrame = min(availableFrames);
    maxFrame = max(availableFrames);
end

if ~isempty(frameRange)
    frameRange = sort(round(frameRange(:).'));
elseif ~isempty(frameStart) && ~isempty(frameWindow)
    frameRange = [round(frameStart), round(frameStart + frameWindow - 1)];
elseif ~isempty(frameWindow)
    frameRange = [minFrame, round(minFrame + frameWindow - 1)];
elseif ~isempty(frameStart)
    frameRange = [round(frameStart), maxFrame];
else
    frameRange = [minFrame, maxFrame];
end

frameMask = MatTracking(:, 4) >= frameRange(1) & MatTracking(:, 4) <= frameRange(2);
points = MatTracking(frameMask, :);
titleText = sprintf('Block %03d localization points | frames %d-%d', ...
    blockIdx, frameRange(1), frameRange(2));

figure;
switch displayMode
    case 'scatter'
        % 散点模式：黑底白色 ROI 横带，半透明红点重叠后会自然加深。
        if isempty(MatTracking)
            displayHeight = 1;
            displayWidth = 1;
            roiZRange = [1, 1];
        else
            roiZRange = [floor(min(MatTracking(:, 2))), ceil(max(MatTracking(:, 2)))];
            displayHeight = max(ceil(max(MatTracking(:, 2))) + 20, roiZRange(2));
            displayWidth = max(ceil(max(MatTracking(:, 3))) + 10, 1);
        end

        backgroundImage = zeros(displayHeight, displayWidth);
        roiZRange(1) = max(1, roiZRange(1));
        roiZRange(2) = min(displayHeight, roiZRange(2));
        backgroundImage(roiZRange(1):roiZRange(2), :) = 1;

        imagesc(backgroundImage);
        colormap(gca, gray);
        hold on;
        if ~isempty(points)
            scatterHandle = scatter(points(:, 3), points(:, 2), pointSize, ...
                'r', 'filled', 'MarkerEdgeColor', 'none');
            scatterHandle.MarkerFaceAlpha = 0.18;
            scatterHandle.MarkerEdgeAlpha = 0.18;
        end
        hold off;
        set(gca, 'YDir', 'reverse');
        axis image;
        xlim([1, displayWidth]);
        ylim([max(1, roiZRange(1) - 20), min(displayHeight, roiZRange(2) + 20)]);

    case 'density'
        % 密度模式：累计定位点并显示为热图，自动聚焦到有点的轴向区域。
        if isempty(points)
            densityImage = zeros(1, 1);
            displayWidth = 1;
            displayHeight = 1;
            zDisplayRange = [1, 1];
        else
            displayWidth = max(ceil(max(MatTracking(:, 3))) + 5, 1);
            displayHeight = max(ceil(max(MatTracking(:, 2))) + 5, 1);

            xBin = max(1, min(displayWidth, round(points(:, 3))));
            zBin = max(1, min(displayHeight, round(points(:, 2))));

            densityMap = accumarray([zBin, xBin], 1, [displayHeight, displayWidth]);

            % 轻微平滑，让连续结构更容易看出来，同时保留点密度差异。
            densityMap = imgaussfilt(densityMap, 0.7);
            densityImage = log10(densityMap + 1);

            % 自动聚焦到当前 block 有定位点的轴向范围，减少上方大片黑区。
            zPadding = 10;
            zDisplayRange = [ ...
                max(1, floor(min(points(:, 2))) - zPadding), ...
                min(displayHeight, ceil(max(points(:, 2))) + zPadding) ...
            ];
        end

        imagesc(1:displayWidth, 1:displayHeight, densityImage);
        set(gca, 'YDir', 'reverse', 'Color', 'k', 'FontSize', 14, 'LineWidth', 1);
        axis image;
        xlim([1, displayWidth]);
        ylim(zDisplayRange);

        colormap(gca, hot(256));
        colorbarHandle = colorbar;
        colorbarHandle.Label.String = 'log10(count + 1)';

        nonzeroValues = sort(densityImage(densityImage > 0));
        if ~isempty(nonzeroValues)
            climIndex = max(1, round(0.995 * numel(nonzeroValues)));
            caxis([0, nonzeroValues(climIndex)]);
        end

        set(gcf, 'Color', 'w', 'Position', [100, 100, 1050, 420]);
end

xlabel('Lateral pixel');
ylabel('Axial pixel');
title(titleText);

fprintf('Block %03d: displayed %d localization points from frames %d-%d.\n', ...
    blockIdx, size(points, 1), frameRange(1), frameRange(2));

end
