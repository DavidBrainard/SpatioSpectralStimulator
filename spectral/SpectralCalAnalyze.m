% SpectralAnalyze
%
% Read and analyze output of SpectralCalCheck
%

% History:
%    11/09/21  dhb  Wrote it.
%    04/19/22  smo  Added contrast residual plot.
%    11/10/22  smo  Added an option to fit and update in format for all
%                   available data at once.

%% Initialize.
clear; close all;

%% Set paramters.
QUICKCHECK = true;
FITALLATONCE = false;
SAVETHEPLOT = true;

%% Fitting happens here.
%
% You can do fit all at once or only single data.
if (FITALLATONCE)
    if (ispref('SpatioSpectralStimulator','SACCData'))
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        testFilenames = dir(fullfile(testFiledir,'CheckCalibration','testImageDataCheck_*'));
    end
    nFits = length(testFilenames);
else
    nFits = 1;
end

% Start a loop here to fit all at once if you want.
for ff = 1:nFits
    %% Load output of SpectralCalCheck
    if (ispref('SpatioSpectralStimulator','SACCData'))
        olderDate = ff-1;
        testFiledir = getpref('SpatioSpectralStimulator','SACCData');
        testFilename = GetMostRecentFileName(fullfile(testFiledir,'CheckCalibration'),...
            'testImageDataCheck_','olderDate',olderDate);
        theCheckData = load(testFilename);
        theData = theCheckData.theData;
    else
        error('Cannot find data file');
    end
    
    % Newly added variables.
    if isfield(theData,{'spatialGaborTargetContrast','targetScreenPrimaryContrast','targetLambda'})
        spatialGaborTargetContrast = theData.spatialGaborTargetContrast;
        targetScreenPrimaryContrast = theData.targetScreenPrimaryContrast;
        targetLambda = theData.targetLambda;
    end
    
    % Get the date of validation.
    [filedir filename ext] = fileparts(testFilename);
    fprintf('\t Current testing file name: (%s) \n', filename);
    numExtract = regexp(filename,'\d+','match');
    yearStr = numExtract{1};
    monthStr = numExtract{2};
    dayStr = numExtract{3};
    hourStr = numExtract{4};
    minuteStr = numExtract{5};
    dateStrVal = sprintf('%s_%s_%s',yearStr,monthStr,dayStr);
    dateStrVal = strrep(dateStrVal,'_','/');
    
    % Get the date of primary measurement.
    if isfield(theCheckData,'primaryFilename')
        [filedir filename ext] = fileparts(theCheckData.primaryFilename);
        numExtract = regexp(filename,'\d+','match');
        dateStrPrimary = sprintf('%s_%s_%s',numExtract{1},numExtract{2},numExtract{3});
        dateStrPrimary = strrep(dateStrPrimary,'_','/');
    end
    
    % We set the target primary contrast as 0.05 till the end of October
    % and it has been upadated to 0.07 after which sounds more sense.
    % However, the performance seems not quite affected by how we set it.
    %
    % For all October data, we used 0.05 for target primary contrasts. And
    % on 11/15 for high contrast test image set, we used 0.10 for target
    % primary contrast. Here, we set it manually, but further data would
    % contain this info in the variable so that we can read out from it.
    if ~exist('targetScreenPrimaryContrast')
        if strcmp(monthStr,'10')
            targetScreenPrimaryContrast = 0.05;
        elseif (strcmp(monthStr,'11') & strcmp(dayStr,'15') & strcmp(hourStr,'16') & strcmp(minuteStr,'13'))
            targetScreenPrimaryContrast = 0.10;
        else
            targetScreenPrimaryContrast = 0.07;
        end
    end
    
    %% Set up some variables that we need
    %
    % Target Spds.  The check here just verifies that we are consistent between
    % current compute code and this measurement.
    targetScreenSpd = theCheckData.targetScreenSpd;
    % targetScreenSpdCompute = theComputeData.screenCalObj.get('P_device');
    % if (any(targetScreenSpd(:) ~= targetScreenSpdCompute(:)))
    %     error('Strange change in target spd');
    % end
    %targetScreenSpd = theComputeData.screenCalObj.get('P_device');
    
    % Set some variables.
    S = theData.S;
    wls = SToWls(S);
    nPrimaries = size(targetScreenSpd,2);
    nChannels = theData.channelCalObjs{1}.get('nDevices');
    channelNInputLevels = size(theData.channelCalObjs{1}.get('gammaInput'),1);
    logicalToPhysical = [0:15];
    nTestPoints = size(theData.ptCldScreenContrastCheckCal,2);
    T_cones = theData.T_cones;
    
    if (~QUICKCHECK)
        %% Load output of SpectralCalCompute.
        %
        % We can remove this now that we are saving testData from SpectralCalCheck,
        % after next round of measurement.
        conditionName = 'LminusMSmooth';
        olderDate = 0;
        if (ispref('SpatioSpectralStimulator','SACCData'))
            testFiledir = getpref('SpatioSpectralStimulator','SACCData');
            testFilename = GetMostRecentFileName(fullfile(testFiledir,'CheckCalibration'),...
                'testImageData_','olderDate',olderDate);
            theComputeData = load(testFilename);
        end
        
        %% Look at measured primaries
        %
        % Make plot comparing what we wanted for primaries versus what we got.
        % What we want is in targetScreenSpd, what we got is in
        % targetScreenSpdMeasured.
        % Plot the spd results.
        primarySpdFig = figure; clf;
        for pp = 1:nPrimaries
            checkTargetPrimaries(:,pp) = SettingsToPrimary(theData.channelCalObjs{pp},theData.screenPrimarySettings(:,pp));
            if (max(abs(checkTargetPrimaries(:,pp)-theData.screenPrimaryPrimaries(:,pp))) > 0.005)
                error('Cannot reconstruct target primary values from target primary settings to quantization tolerance');
            end
            checkTargetSettings(:,pp) = PrimaryToSettings(theData.channelCalObjs{pp},theData.screenPrimaryPrimaries(:,pp));
            if (any(theData.screenPrimarySettings(:,pp) ~= checkTargetSettings(:,pp)))
                error('Cannot reconstruct target primary settings from calibraion and primary values');
            end
            checktargetScreenSpd(:,pp) = PrimaryToSpd(theData.channelCalObjs{pp},checkTargetPrimaries(:,pp));
            if (any(targetScreenSpd(:,pp) ~= checktargetScreenSpd(:,pp)))
                error('Cannot reconstruct target primary spd from calibraion and primary values');
            end
            
            subplot(nPrimaries,1,pp); hold on;
            plot(wls,targetScreenSpd(:,pp),'k','LineWidth',4)
            plot(wls,theCheckData.targetScreenSpdMeasured(:,pp),'r','LineWidth',3);
            plot(wls,checktargetScreenSpd(:,pp),'g','LineWidth',1);
            legend('Target','Measured','Target Check');
            xlabel('Wavelength (nm)');
            ylabel('Spectral power distribution');
            title('Comparison of raw measured and desired spds');
        end
        
        %% Use regression to express measured primaries in terms of channel calibration
        primaryValueFig = figure; clf;
        figureSize = 1000;
        figurePosition = [1200 300 figureSize figureSize/3];
        set(gcf,'position',figurePosition);
        for pp = 1:nPrimaries
            P_device{pp} = theData.channelCalObjs{pp}.get('P_device');
            P_ambient{pp} = theData.channelCalObjs{pp}.get('P_ambient');
            if (any(P_ambient{pp}(:) ~= 0))
                error('Need to handle non-zero ambient');
            end
            regressTargetPrimaries(:,pp) = P_device{pp}\theCheckData.targetScreenSpdMeasured(:,pp);
            
            figure(primaryValueFig);
            subplot(1,nPrimaries,pp); hold on;
            plot(checkTargetPrimaries(:,pp),regressTargetPrimaries(:,pp),'ro','MarkerSize',12,'MarkerFaceColor','r');
            xlim([0,1]); ylim([0,1]); axis('square');
            plot([0 1],[0 1],'k');
            title('Regression versus desired channel values')
            xlabel('Desired Primary'); ylabel('Regression on Measured Primary');
            
            figure(primarySpdFig);
            subplot(nPrimaries,1,pp); hold on;
            plot(wls,P_device{pp}*regressTargetPrimaries(:,pp),'b','LineWidth',1);
            legend('Target','Measured','Target Check','Regress');
        end
        
        %% Make plot of measured versus desired spds.
        %
        % We should be able to obtain these as a linear combination
        % of the primaries.
        plotRegress = false;
        figure; clf;
        figureSize = 1000;
        figurePosition = [1200 300 figureSize figureSize];
        set(gcf,'position',figurePosition);
        for tt = 1:nTestPoints
            % Find linear combination of measured primaries to produce test
            regressPointCloudPrimaries(:,tt) = theCheckData.targetScreenSpdMeasured\theCheckData.ptCldScreenSpdMeasuredCheckCal(:,tt);
            regressPointCloudSpd(:,tt) = theCheckData.targetScreenSpdMeasured*regressPointCloudPrimaries(:,tt);
            
            subplot(round(nTestPoints/2),2,tt); hold on;
            plot(wls,theCheckData.ptCldScreenSpdCheckCal(:,tt),'k-','LineWidth',4)  % Target spectra
            plot(wls,theCheckData.ptCldScreenSpdMeasuredCheckCal(:,tt),'r-','LineWidth',2); % Measured spectra
            xlabel('Wavelength (nm)');
            ylabel('Spectral power distribution');
            title(sprintf('Test %d raw',tt),'fontsize',16);
            if (plotRegress)
                plot(wls,regressPointCloudSpd(:,tt),'g-','LineWidth',1);              % Regression fit
                legend('Target','Measured','Regress');
            else
                legend('Target','Measured');
            end
        end
        
        % Another way of comparing measured and nominal spectra
        ptCldSpdScatter = theCheckData.ptCldScreenSpdMeasuredCheckCal;
        figure; clf;
        figureSize = 1000;
        figurePosition = [1200 300 figureSize figureSize];
        set(gcf,'position',figurePosition);
        spdLim = 2.5e-3;
        for tt = 1:nTestPoints
            scaleFactor(tt) = theCheckData.ptCldScreenSpdCheckCal(:,tt)\theCheckData.ptCldScreenSpdMeasuredCheckCal(:,tt);
            subplot(round(nTestPoints/2),2,tt); hold on;
            plot(theCheckData.ptCldScreenSpdCheckCal(:,tt),theCheckData.ptCldScreenSpdMeasuredCheckCal(:,tt),'r+');
            xlabel('Nominal Spd');
            ylabel('Measured Spd');
            xlim([0 spdLim]); ylim([0 spdLim]);
            axis('square');
            title(sprintf('Test %d raw, factor %0.3f',tt,scaleFactor(tt)),'fontsize',12)
        end
        
        % Here we will compare the primaries we get from regression above to those
        % we wanted, and see how we do.
        figure; clf;
        figureSize = 1000;
        figurePosition = [1200 300 figureSize figureSize/3];
        set(gcf,'position',figurePosition);
        for pp = 1:nPrimaries
            subplot(1,nPrimaries,pp); hold on;
            plot(theCheckData.ptCldScreenPrimariesCheckCal(pp,:),regressPointCloudPrimaries(pp,:),'ro','MarkerSize',12,'MarkerFaceColor','r');
            plot([0 1],[0 1],'k');
            xlim([0,1]); ylim([0,1]); axis('square');
            title('Regression versus desired primary values')
            xlabel('Desired Screen Primary'); ylabel('Regression on Measured Screen Primary');
        end
        
        %% Explicitly compute some contrasts
        analyzePointCloudExcitationsMeasured = T_cones * theCheckData.ptCldScreenSpdMeasuredCheckCal;
        analyzePointCloudBgExcitationsMeasured = analyzePointCloudExcitationsMeasured(:,1);
        analyzePointCloudContrastMeasured = ExcitationToContrast(analyzePointCloudExcitationsMeasured,analyzePointCloudBgExcitationsMeasured);
        if (any(analyzePointCloudContrastMeasured(:) ~= theCheckData.ptCldScreenContrastMeasuredCheckCal(:)))
            error('Cannot get same measured contrasts in two different places');
        end
        
        analyzePointCloudExcitationsNominal = T_cones * theCheckData.ptCldScreenSpdCheckCal;
        analyzePointCloudBgExcitationsNominal = analyzePointCloudExcitationsNominal(:,1);
        analyzePointCloudContrastNominal = ExcitationToContrast(analyzePointCloudExcitationsNominal,analyzePointCloudBgExcitationsNominal);
        if (any(analyzePointCloudContrastNominal(:) ~= theCheckData.ptCldContrastNominal(:)))
            error('Cannot get same nominalcontrasts in two different places');
        end
    end
    
    %% Plot measured versus desired contrasts
    contrastFig = figure; hold on;
    figureSize = 1100;
    figurePosition = [1200 300 figureSize figureSize/3];
    set(gcf,'position',figurePosition);
    
    if (theData.spatialGaborTargetContrast > 0.07)
        axisLim = 0.10;
    else
        axisLim = 0.07;
    end
    
    theColors = ['r' 'g' 'b'];
    for pp = 1:nPrimaries
        subplot(1,nPrimaries,pp); hold on;
        markerSizeFilled = 19;
        markerSizeUnfilled = 22;
        plot(theCheckData.desiredContrastCheckCal(pp,:),theCheckData.ptCldScreenContrastMeasuredCheckCal(pp,:),[theColors(pp) 'o'],'MarkerSize',markerSizeFilled,'MarkerFaceColor',theColors(pp));
        plot(theCheckData.desiredContrastCheckCal(pp,:),theCheckData.ptCldContrastNominal(pp,:), [theColors(pp) 'o'],'MarkerSize',markerSizeUnfilled);
        plot(theCheckData.desiredContrastCheckCal(pp,1),theCheckData.ptCldScreenContrastMeasuredCheckCal(pp,1),'ko','MarkerSize',markerSizeFilled,'MarkerFaceColor','k');
        plot(theCheckData.desiredContrastCheckCal(pp,1),theCheckData.ptCldContrastNominal(pp,1), 'ko','MarkerSize',markerSizeUnfilled);
        
        plot([-1 1],[-1 1],'k');
        xlim([-axisLim axisLim]);
        ylim([-axisLim axisLim]);
        axis('square');
        fontSize = 15;
        xlabel('Desired cone contrast','fontsize',fontSize);
        ylabel('Measured cone contrast','fontsize',fontSize);
        legend({'Measured','Nominal'},'location','southeast','fontsize',fontSize);
        title(sprintf('Cone class %d',pp),'fontsize',fontSize);
    end
    
    % Set the primary measurement date and validation date right. For most
    % cases, validation was performed on the same day when the primary
    % measurement was made except some cases. Here we added the cases so
    % that it shows the right dates for both.
    if ~exist('dateStrPrimary')
        if (strcmp(monthStr,'11') & strcmp(dayStr,'06'))
            if (strcmp(hourStr,'17') & strcmp(minuteStr,'13'))
                dateStrPrimary = '2022/10/31';
                targetScreenPrimaryContrast = 0.05;
            elseif (strcmp(hourStr,'17') & strcmp(minuteStr,'06'))
                dateStrPrimary = '2022/11/03';
            elseif (strcmp(hourStr,'16') & strcmp(minuteStr,'57'))
                dateStrPrimary = '2022/10/26';
                targetScreenPrimaryContrast = 0.05;
            end
        else
            dateStrPrimary = dateStrVal;
        end
    end
    
    % Add some texts to the plot.
    main = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
    text(0.08,0.97,sprintf('* Date of primary measurement: ( %s )',dateStrPrimary),'fontsize',15,'Parent',main);
    text(0.08,0.92,sprintf('* Date of validation: ( %s )',dateStrVal),'fontsize',15,'Parent',main);
    text(0.40,0.95,sprintf('* Target Primary Contrast: ( %.2f )',targetScreenPrimaryContrast),'fontsize',15,'Parent',main);
    text(0.68,0.95,sprintf('* Target Image Contrast: ( %.2f )',theData.spatialGaborTargetContrast),'fontsize',15,'Parent',main);
    
    if (SAVETHEPLOT)
        if (ispref('SpatioSpectralStimulator','SACCAnalysis'))
            testFiledir = fullfile(getpref('SpatioSpectralStimulator','SACCAnalysis'),'CheckCalibration');
            
            % Back to underbar to save it on the file name.
            dateStrVal = strrep(dateStrVal,'/','_');
            dateStrPrimary = strrep(dateStrPrimary,'/','_');
            
            % Save the plot.
            testFilename = fullfile(testFiledir,sprintf('testImageDataCheck_P(%s)_V(%s)_PC(%.2f)_IC(%.2f)',...
                dateStrPrimary,dateStrVal,targetScreenPrimaryContrast,theData.spatialGaborTargetContrast));
            testFileFormat = '.tiff';
            saveas(gcf,append(testFilename,testFileFormat));
            fprintf('\t Plot has been saved successfully! \n');
        end
    end
    
    %% Plot the residual between desired and measured contrasts.
    if (~QUICKCHECK)
        contrastResidualFig = figure; hold on;
        figureSize = 1000;
        figurePosition = [1200 300 figureSize figureSize/3];
        set(gcf,'position',figurePosition);
        
        xAxisLim = 0.03;
        yAxisLim = 0.008;
        theColors = ['r' 'g' 'b'];
        contrastResidual = theCheckData.ptCldScreenContrastMeasuredCheckCal -theCheckData.desiredContrastCheckCal;
        
        for pp = 1:nPrimaries
            subplot(1,nPrimaries,pp); hold on;
            
            % All contrast test points.
            desiredContrast = theCheckData.desiredContrastCheckCal;
            
            % As the desired S cone contrasts are all zeros, so we set this field
            % with desired L cone contrast.
            desiredContrast(3,:) = desiredContrast(1,:);
            
            % Plot the contrast residual here.
            plot(desiredContrast(pp,:),contrastResidual(pp,:),[theColors(pp) 'o'],'MarkerSize',14,'MarkerFaceColor',theColors(pp));
            plot(desiredContrast(pp,:),zeros(1,size(theCheckData.desiredContrastCheckCal,2)),[theColors(pp) 'o'],'MarkerSize',18);
            
            % Plot the zero contrast for the reference here.
            plot(desiredContrast(pp,1),contrastResidual(pp,1),'ko','MarkerSize',14,'MarkerFaceColor','k');
            plot(desiredContrast(pp,1),zeros(1,1),'ko','MarkerSize',18);
            
            % Draw the reference horizontal line.
            plot([-1 1],[0 0],'k');
            
            xlim([-xAxisLim xAxisLim]);
            ylim([-yAxisLim yAxisLim]);
            axis('square');
            xlabel('Desired contrast');
            ylabel('Contrast Residual (Measured - Desired)');
            legend({'Measured','Nominal'},'location','southeast');
            title(sprintf('Cone class %d',pp));
        end
    end
end
