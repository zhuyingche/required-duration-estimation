%----function to combine nmfEpochResult of different states into one struct (nmfEpochResults)--%
%20250616: 1. Current order: InterictalDataOnly, RmNan,
% RmFewValidSampleEpOption, RmAbnEpOption; 
% previous order: RmFewValidSampleEpOption,RmAbnEpOption,
% RmNan,InterictalDataOnly;
% 2. RmAbnEpOption: change rmLow criteria
function combineDiffStateIntoOne(curPt,varNames, filePathResult, hfoversion, winMin, RmAbnEpOption,RmFewValidSampleEpOption,normalizeMethod,InterictalDataOnly)

resultPath = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'], 'nmf', RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod,InterictalDataOnly);%before 2025 & 20251010

    for i = 1:numel(varNames)
        varName = varNames{i};
        % resultPath = fullfile(filePathResult, hfoversion,
        % [num2str(winMin) 'min'], 'nmf', InterictalDataOnly,
        % RmFewValidSampleEpOption,RmAbnEpOption,normalizeMethod);%20250616
        resultFile = fullfile(resultPath, strcat(curPt, '.nmfEpochResult.', varName, '.mat'));
        if isfile(resultFile)
            load(resultFile);
            nmfEpochResults.(varName) = nmfEpochResult;
            clear nmfEpochResult
        end
    end

     % Save the data
     nameF = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'],'nmf',RmAbnEpOption,RmFewValidSampleEpOption,normalizeMethod,InterictalDataOnly);%before 2025 & 20251010
        % nameF = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'],'nmf',InterictalDataOnly, RmFewValidSampleEpOption,RmAbnEpOption,normalizeMethod);%20250616
        if ~exist(nameF, 'dir')
            mkdir(nameF);
        end
        filename = fullfile(nameF, strcat(curPt, '.nmfEpochResult.mat'));
        % save(filename, 'nmfEpochResults');

       if exist("fileStartTimeMin") && exist("fileStopTimeMin")
            save(filename, 'nmfEpochResults','fileStartTimeMin','fileStopTimeMin');
        else
            save(filename, 'nmfEpochResults');
       end


   % ---- 20251010: Delete individual files after successful save ----
    for i = 1:numel(varNames)
        resultPath = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'], ...
            'nmf', RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly);
        resultFile = fullfile(resultPath, strcat(curPt, '.nmfEpochResult.', varNames{i}, '.mat'));

        if isfile(resultFile)
            try
                delete(resultFile);
            catch ME
                warning('Could not delete file: %s\nReason: %s', resultFile, ME.message);
            end
        end
    end

end