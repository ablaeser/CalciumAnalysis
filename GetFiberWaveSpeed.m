function [fiber, pooledFit] = GetFiberWaveSpeed(expt, ROI, fiber, periBout, varargin) % , waveType
checkBout = @(x)(isstruct(x) || iscell(x));
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'ROI', @isstruct )
addRequired( IP, 'fiber', @isstruct )
addRequired( IP, 'periBout', checkBout ) %  @isstruct
addOptional( IP, 'csdWave', [], @iscell) %
addParameter( IP, 'lp', 2, @isnumeric) % cutoff frequency (Hz) for Gaussian LP filter
addParameter( IP, 'minLength', 200, @isnumeric)
addParameter( IP, 'minEffect', 1, @isnumeric)
addParameter( IP, 'minR2', 0.8, @isnumeric )
addParameter( IP, 'maxTau', 2, @isnumeric )
addParameter( IP, 'show', false, @islogical)
parse( IP, expt, ROI, fiber, periBout, varargin{:} );  % , waveType
csdWave = IP.Results.csdWave; csdWave = csdWave{1};
lpFreq = IP.Results.lp;
minLength = IP.Results.minLength;
minEffect  = IP.Results.minEffect;
minR2 = IP.Results.minR2;
maxTau = IP.Results.maxTau;
show = IP.Results.show;

if lpFreq > expt.scanRate, lpFreq = expt.scanRate/(2*pi); end
gaussSigma = 1/(2*pi*lpFreq);
gaussFilt = MakeGaussFilt( 4, 0, gaussSigma, expt.scanRate, false ); %  gaussFilt = MakeGaussFilt( 4, 0, 1/expt.scanRate, expt.scanRate, true );
Nfiber = numel(fiber);
fiberColor = distinguishable_colors(Nfiber);
[fiberLengthSort, fSort] = sort([fiber.length], 'descend');

% Set up sigmoidal step function fit
modelStep = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4)./(1+exp((x-b(6))/b(5))) + b(7)); %
stepInit(1) = 1; % onset amplitude
stepInit(2) = 0.1; % onset time constant
stepInit(3) = 0; % Onset time
stepInit(4) = 1; % offset amplitude
stepInit(5) = 0.2; % offset time constant
stepInit(6) = 3; % Offset time
stepInit(7) = 0; % Constant term
lsqfitOpts = struct('Display','off'); % fitOpts = statset('Display','off', 'MaxIter',1000);

% Set up sigmoidal fit
%{
modelSigmoid = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4) ); % beta = [factor, growth rate, threshold, constant response offset]
sigmoidInit(1) = 1; % Saturation response
sigmoidInit(2) = 0.1; % time constant
sigmoidInit(3) = 0; % Onset time
sigmoidInit(4) = 0; % Constant offset response
%}
opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
if show, close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]); end
if ~isempty(csdWave)
    waveType = 'CSD';
    onsetPool = zeros(0,2);
    for b = 1
        for f = fSort(fiberLengthSort > minLength) % 1:Nfiber
            % Get peri-CSD fiber fluor signals
            %imagesc( periBout.fluor{b}(csdWave.scans,fiber(f).ROI )' ); impixelinfo;
            fiber(f).(waveType).Twindow = csdWave.T;
            for r = 1:fiber(f).Nroi
                % Which bin does this ROI belong to?
                tempCent = ROI(fiber(f).ROI(r)).cent;
                xBinTemp = find( tempCent(1) > csdWave.xBinLim, 1, 'last');
                if xBinTemp > csdWave.NbinX, xBinTemp = csdWave.NbinX; end
                if isempty(xBinTemp), xBinTemp = 1; end
                yBinTemp = find( tempCent(2) > csdWave.yBinLim, 1, 'last');
                if isempty(yBinTemp), yBinTemp = 1; end
                if yBinTemp > csdWave.NbinY, yBinTemp = csdWave.NbinY; end
                tempBinCent = squeeze(csdWave.binCent(xBinTemp, yBinTemp,1:2))'/expt.umPerPixel;
                fiber(f).(waveType).binOnset(b,r) = csdWave.binOnsetTime(xBinTemp,yBinTemp);
                
                % Normalize the ROI's trace
                tempTrace = periBout.fluor{b}(csdWave.scans,fiber(f).ROI(r)); % :
                %check for missing scans
                nanScan = isnan(tempTrace);
                Nnan = sum(nanScan);
                fprintf('\n%i missing values found. %i good values remain', Nnan, sum(~nanScan) )
                %{
                if  Nnan > 0
                    % drop leading nan if possible
                    tempCC = bwconncomp( nanScan );
                    dropScan = tempCC.PixelIdxList{1};
                    if ismember(1, dropScan)
                        dropScan = dropScan(1:csdWave.scans(1)-1);
                        tempTrace(tempCC.PixelIdxList{1}) = [];
                        fprintf('\n%i missing values found. Removed %i leading missing values, %i values remain', Nnan, numel(tempCC.PixelIdxList{1}), numel(tempTrace) )
                    end
                end
                %}
                try
                    tempTrace = filtfilt(gaussFilt, 1, tempTrace); % plot(periBout.fluor{b}(:,r)); hold on; plot(tempTrace); hold on;
                catch
                    fprintf('\n[b,f,r] = [%i, %i, %i]: Filtering failed', b, f, r);
                end
                %tempTrace = tempTrace(csdWave.scans,:);
                tempEffect = mean(tempTrace(csdWave.postScan), 'omitnan') - median(tempTrace(csdWave.preScan), 'omitnan'); % max(tempTrace(csdWave.postScan),[], 'omitnan')
                tempTrace = (tempTrace - prctile(tempTrace,10)); % median(tempTrace(csdWave.preScan), 'omitnan')
                tempTrace = tempTrace./max(tempTrace(csdWave.postScan));
                fiber(f).(waveType).trace(:,r) = tempTrace;
                
                % Fit peri-onset signal to step
                if ~isnan(csdWave.binOnsetTime(xBinTemp,yBinTemp))
                    stepInit(2) = csdWave.binTimeConst(xBinTemp,yBinTemp); % onset time constant
                    stepInit(3) = csdWave.binOnsetTime(xBinTemp,yBinTemp); % Onset time
                    stepInit(5) = csdWave.binTimeConst(xBinTemp,yBinTemp); % offset time constant
                    stepInit(6) = 3; % Offset time
                else
                    stepInit(2) = median(csdWave.binTimeConst(:), 'omitnan'); % onset time constant
                    stepInit(3) = 0; % Onset time
                    stepInit(5) = median(csdWave.binTimeConst(:), 'omitnan'); % offset time constant
                    stepInit(6) = 3; % Offset time
                end
                fitResult = struct('R2',NaN, 'Tonset',NaN, 'tauOn',NaN, 'Toffset',NaN, 'tauOff',NaN, 'pTonset',NaN, 'pTauOn',NaN, 'pToffset',NaN, 'pTauOff',NaN, 'prediction',nan(length(tempTrace),1)); % struct('R2',NaN, 'Tonset',NaN, 'pTonset',NaN, 'tau',NaN, 'pTau',NaN);
                % Set upper and lower bounds on the parameters: important to allow for possibility of long delay in onset time
                %[Onset amplitude, onset time constant, Onset time, offset amplitude, offset time constant, Offset time, Constant term]
                stepLower = [0, 0,      -10, 0,      0,  0,                                -1];
                stepUpper = [2, maxTau,  10, maxTau, 2,  fiber(f).(waveType).Twindow(end),  1];
                if ~all(isnan(tempTrace)) && tempEffect > 0
                    %fiber(f).(waveType).traceFit{r} = fitnlm( fiber(f).(waveType).Twindow, tempTrace, modelSigmoid, sigmoidInit, 'Options',fitOpts );
                    [stepFit, resnorm] = lsqcurvefit(modelStep, stepInit, fiber(f).(waveType).Twindow, tempTrace, stepLower, stepUpper, lsqfitOpts ); % , resid, exitFlag
                    fitResult.R2 = 1-resnorm/sum((tempTrace - mean(tempTrace)).^2); %R2; %fiber(f).(waveType).traceFit{k,roi}.Rsquared.Adjusted; % goodness of fit
                    fitResult.Tonset = stepFit(3); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(3);
                    fitResult.tauOn = stepFit(2); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(2);
                    fitResult.Toffset = stepFit(6); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(6);
                    fitResult.tauOff = stepFit(5); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(5);
                    fitResult.prediction = modelStep(stepFit, fiber(f).(waveType).Twindow);
                else
                    fprintf('\nSkipped fitting - no data available');
                end
                fiber(f).(waveType).R2(b,r) = fitResult.R2;
                fiber(f).(waveType).effect(b,r) = tempEffect;
                if tempEffect > minEffect && fitResult.R2 > minR2 %&& fitResult.Toffset >= fitResult.Tonset
                    fiber(f).(waveType).goodFit(b,r) = true;
                else
                    fiber(f).(waveType).goodFit(b,r) = false;
                end
                if fiber(f).(waveType).goodFit(b,r) %tempEffect > minEffect && fitResult.R2 > minR2 && fitResult.tau < maxTau && fitResult.tau >= 0
                    fiber(f).(waveType).tau(b,r) = fitResult.tauOn;
                    fiber(f).(waveType).Tonset(b,r) = fitResult.Tonset;
                    fiber(f).(waveType).onsetDelay(b,r) = fiber(f).(waveType).Tonset(b,r) - csdWave.binOnsetTime(xBinTemp,yBinTemp);
                else
                    fiber(f).(waveType).tau(b,r) = NaN;
                    fiber(f).(waveType).Tonset(b,r) = NaN;
                    fiber(f).(waveType).onsetDelay(b,r) = NaN;
                end
                
                if show
                    subtightplot(2,2,1,opt{:}); cla;
                    imshow( label2rgb(fiber(f).labelFoot), []);
                    
                    subtightplot(2,2,3,opt{:}); cla;
                    %imagesc(csdWave.binOnsetTime'); colorbar;
                    imshow(max(ROI(fiber(f).ROI(r)).labelVol,[],3), []); hold on;
                    plot( tempCent(1), tempCent(2), '.', 'color',fiberColor(f,:))
                    for x = 1:csdWave.NbinX+1,  line(csdWave.xBinLim(x)*[1,1], [1,expt.Nrow], 'color','w' );  end
                    for y = 1:csdWave.NbinY+1,  line([1,expt.Ncol], csdWave.yBinLim(y)*[1,1], 'color','w');   end
                    for x = 1:csdWave.NbinX
                        for y = 1:csdWave.NbinY
                            plot( csdWave.binCent(x,y,1)/expt.umPerPixel, csdWave.binCent(x,y,2)/expt.umPerPixel, 'w.');
                        end
                    end
                    plot( tempBinCent(1), tempBinCent(2), '.', 'color',fiberColor(f,:))
                    
                    subtightplot(2,2,2,opt{:}); cla;
                    plot( csdWave.T, squeeze(csdWave.normTrace(xBinTemp,yBinTemp,:)), 'k' ); hold on;
                    plot( csdWave.T, tempTrace );
                    plot( csdWave.T, fitResult.prediction, 'r', 'lineStyle','--' ); %try  end % predict(fiber(f).(waveType).traceFit{r}, fiber(f).(waveType).Twindow)
                    legend('CSD', 'ROI', 'Fit', 'Location','northWest');
                    axis square;
                    ylim([-Inf,Inf]); % [-0.02, 1.02]
                    xlabel('Peri-onset time (s)'); ylabel('Normalized Fluor');
                    title(sprintf('[f,r] = [%i, %i] effect = %2.2f, R^2 = %2.2f, tau = %2.2f, onset = %2.2f. Difference = %2.2f s', ...
                        f, r, tempEffect, fitResult.R2, fitResult.tauOn, fitResult.Tonset, fitResult.Tonset-csdWave.binOnsetTime(xBinTemp,yBinTemp) ) );
                    
                    subtightplot(2,2,4,opt{:});
                    line([0,10], [0,10], 'color','r', 'lineStyle','--'); hold on;
                    plot( fiber(f).(waveType).binOnset(b,r), fiber(f).(waveType).Tonset(b,r), '.', 'color',fiberColor(f,:) ); hold on;
                    xlabel('Neuropil Onset Time (s)'); ylabel('ROI Onset Time (s)');
                    axis square;
                    pause;
                end
            end
            % Linear fit - neuropil vs ROI onset
            fiber(f).(waveType).linFit = fitlm( fiber(f).(waveType).binOnset', fiber(f).(waveType).Tonset' ); % fit individual fibers
            if show, plot( fiber(f).(waveType).binOnset', predict(fiber(f).(waveType).linFit, fiber(f).(waveType).binOnset'), 'lineStyle','--', 'color',fiberColor(f,:) ); end
            onsetPool = [onsetPool; [fiber(f).(waveType).binOnset', fiber(f).(waveType).Tonset']];
            
            % Fit pairwise separation vs delay to a linear function
            tempROI = find(~all(isnan(fiber(f).(waveType).Tonset),1));
            fiber(f).(waveType).sepDelay = nan(0, 2);
            if numel(tempROI) > 1
                tempPair = nchoosek(tempROI, 2);
                for p = 1:size(tempPair,1)
                    pairSep = fiber(f).footSep(tempPair(p,1), tempPair(p,2));
                    pairDelays = diff(fiber(f).(waveType).Tonset(:,[tempPair(p,1), tempPair(p,2)]),1,2);
                    fiber(f).(waveType).sepDelay = [fiber(f).(waveType).sepDelay; [repmat(pairSep,size(pairDelays,1),1), pairDelays]];
                end
            end
            if sum(~isnan(fiber(f).(waveType).sepDelay(:,2))) > 1
                fiber(f).(waveType).sepDelayFit = fitlm(fiber(f).(waveType).sepDelay(:,1), fiber(f).(waveType).sepDelay(:,2), 'RobustOpts',true);
                if fiber(f).(waveType).sepDelayFit.Coefficients.pValue(1) < 0.05
                    fiber(f).(waveType).sepDelayConst = fiber(f).(waveType).sepDelayFit.Coefficients.Estimate(1);
                else
                    fiber(f).(waveType).sepDelayConst = 0;
                end
                if fiber(f).(waveType).sepDelayFit.Coefficients.pValue(2) < 0.05
                    fiber(f).(waveType).sepDelaySlope = fiber(f).(waveType).sepDelayFit.Coefficients.Estimate(2);
                else
                    fiber(f).(waveType).sepDelaySlope = 0;
                end
                fiber(f).(waveType).sepDelayPred = polyval([fiber(f).(waveType).sepDelaySlope; fiber(f).(waveType).sepDelayConst], fiber(f).(waveType).sepDelay(:,1));
            else
                fiber(f).(waveType).sepDelayFit = [];
                fiber(f).(waveType).sepDelayConst = NaN;
                fiber(f).(waveType).sepDelaySlope = NaN;
                fiber(f).(waveType).sepDelayPred = [];
            end
            if show && ~isempty(fiber(f).(waveType).sepDelayFit)
                clf;
                plot( fiber(f).(waveType).sepDelay(:,1), fiber(f).(waveType).sepDelay(:,2), '.'); hold on;
                plot( fiber(f).(waveType).sepDelay(:,1), fiber(f).(waveType).sepDelayPred, 'r--' )
                xlabel('Pairwise Separation (um)'); ylabel('Pairwise Lag (s)');
                title(sprintf('CSD Onset: f = %i: ', f));
                axis square;
                pause;
            end
        end
        % {
        onsetPool = onsetPool(~isnan(sum(onsetPool,2)),:);
        if size(onsetPool,1) >= 2 && show
            pooledFit = fitlm( onsetPool(:,1), onsetPool(:,2) ); % fit pooled fibers  , 'RobustOpts',true
            %line([0,10], [0,10], 'color','r', 'lineStyle','--'); hold on;
            plot(onsetPool(:,1), onsetPool(:,2), '.'); hold on; axis square;
            plot(onsetPool(:,1), predict(pooledFit, onsetPool(:,1)), 'r', 'lineStyle','--' );
        end
        %}
    end
else
    waveType = 'spont';
    for f = fSort(fiberLengthSort > minLength) %  1:Nfiber
        fiber(f).(waveType).runBout = nan(0,2); % indices (run,b) for each analyzed bout
        fiber(f).(waveType).Twindow = cell(0,0); fiber(f).(waveType).trace = cell(0,0); fiber(f).(waveType).traceFit = cell(0,0); fiber(f).(waveType).goodFit = false(0,fiber(f).Nroi);
        fiber(f).(waveType).effect = nan(0, fiber(f).Nroi); fiber(f).(waveType).R2 = nan(0, fiber(f).Nroi);
        fiber(f).(waveType).tauOn = nan(0, fiber(f).Nroi); fiber(f).(waveType).Tonset = nan(0, fiber(f).Nroi);
        k = 0;
        %periBout(run).stat.fluor.effect(:,fiber(f).ROI,1) > 1
        for run = 1:numel(periBout)
            effectCheck = periBout(run).stat.fluor.effect(:,fiber(f).ROI,1) > minEffect;
            bEffect = find( sum(effectCheck,2) > 1 )'; % must have at least 2 ROI showing effect
            [effectSort, bSort] = sort( median( periBout(run).stat.fluor.effect(:,fiber(f).ROI,1), 2 ), 'descend', 'MissingPlacement','last' );
            for b = intersect(bSort', bEffect, 'stable') % %1:periBout(run).Nbout % (effectSort > minEffect)' %
                if max( fiber(f).footSep(effectCheck(b,:), find(effectCheck(b,:),1,'first')) ) > minLength % make sure the good ROIs are sufficiently far apart to be useful
                    k = k+1;
                    fiber(f).(waveType).runBout(k,:) = [run, b];
                    Twindow = periBout(run).T{b} - periBout(run).Tstart(b); % (windowScans)
                    fiber(f).(waveType).Twindow{k} = Twindow;
                    fiber(f).(waveType).Tonset(k,1:fiber(f).Nroi) = NaN;
                    for roi = flip(find(effectCheck(b,:))) %flip(1:fiber(f).Nroi)
                        fiber(f).(waveType).effect(k,roi) = periBout(run).stat.fluor.effect(b,fiber(f).ROI(roi),1);
                        tempTrace = periBout(run).fluor{b}(:,fiber(f).ROI(roi)); % filtfilt(gaussFilt, 1, periBout.fluor{b}(:,fiber(f).ROI(roi)) ); %
                        Nnan = sum(isnan(tempTrace));
                        fitResult = struct('R2',NaN, 'Tonset',NaN, 'tauOn',NaN, 'Toffset',NaN, 'tauOff',NaN, 'pTonset',NaN, 'pTauOn',NaN, 'pToffset',NaN, 'pTauOff',NaN, 'prediction',nan(length(tempTrace),1));
                        if Nnan == 0
                            % Normalize the trace
                            tempTrace = (tempTrace - median( tempTrace([periBout(run).preScan{b}; periBout(run).postScan{b}]), 'omitnan')); % ,
                            tempTrace = tempTrace/max(tempTrace, [], 1, 'omitnan');
                            fiber(f).(waveType).trace{k,roi} = tempTrace;
                            %plot(tempTrace)  %plot(periBout.fluor{b}(:,roi)); hold on; plot(tempFiltTrace);
                            
                            % Fit peri-onset signal to step function
                            stepInit(6) = periBout(run).dur(b); % Offset time
                            try
                                % Set upper and lower bounds on the parameters: important to allow for possibility of long delay in onset time
                                %[Onset amplitude, onset time constant, Onset time, offset amplitude, offset time constant, Offset time, Constant term]
                                stepLower = [0, 0,      -10, 0,      0,  0,            -1];
                                stepUpper = [2, maxTau,  10, maxTau, 2,  Twindow(end),  1];
                                [stepFit, resnorm] = lsqcurvefit(modelStep, stepInit, Twindow, tempTrace, stepLower, stepUpper, lsqfitOpts ); % , resid, exitFlag
                                %fiber(f).(waveType).traceFit{k,roi} = fitnlm( Twindow, tempTrace, modelStep, stepInit, 'Options',fitOpts ); % modelSigmoid
                                fitResult.R2 = 1-resnorm/sum((tempTrace - mean(tempTrace)).^2); %fiber(f).(waveType).traceFit{k,roi}.Rsquared.Adjusted; % goodness of fit
                                fitResult.Tonset = stepFit(3); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(3);
                                fitResult.tauOn = stepFit(2); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(2);
                                fitResult.Toffset = stepFit(6); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(6);
                                fitResult.tauOff = stepFit(5); %fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(5);
                                fitResult.prediction = modelStep(stepFit, fiber(f).(waveType).Twindow{k});
                            catch
                                fprintf('Skipped fitting - no data available');
                            end
                            
                            % Save fit results to fiber struct
                            fiber(f).(waveType).R2(k,roi) = fitResult.R2;
                            if fitResult.R2 > minR2 && fiber(f).(waveType).effect(k,roi) > minEffect %fitResult.tauOn <= maxTau && fitResult.tauOn >= 0 %  &&
                                fiber(f).(waveType).goodFit(k,roi) = true;
                            else
                                fiber(f).(waveType).goodFit(k,roi) = false;
                            end
                            if fiber(f).(waveType).goodFit(k,roi) %fitResult.R2 > minR2 && fitResult.tauOn < maxTau && fitResult.tauOn >= 0 && fiber(f).(waveType).effect(k,roi) > minEffect
                                fiber(f).(waveType).tauOn(k,roi) = fitResult.tauOn;
                                fiber(f).(waveType).Tonset(k,roi) = fitResult.Tonset;
                            else
                                fiber(f).(waveType).tauOn(k,roi) = NaN;
                                fiber(f).(waveType).Tonset(k,roi) = NaN;
                            end
                            if show %&& fitResult.R2 > minR2 && fitResult.tauOn < maxTau && fitResult.tauOn >= 0 && fiber(f).(waveType).effect(k,roi) > minEffect
                                subtightplot(3,1,1,opt{:}); cla;
                                imagesc( periBout(run).fluor{b}(:,fiber(f).ROI )' ); hold on;
                                plot(interp1(Twindow, 1:numel(Twindow), fiber(f).(waveType).Tonset(k,roi) ), roi, 'ko')
                                impixelinfo; 
                                
                                set(gca, 'Ytick',1:fiber(f).Nroi);
                                title(sprintf('[fiber, run, bout] = [%i, %i, %i]', f, run, b) );
                                
                                subtightplot(3,1,2,opt{:}); cla;
                                plot( Twindow, tempTrace ); hold on;
                                plot( Twindow, fitResult.prediction, 'r', 'lineStyle','--' );
                                %plot( periBout(run).T{b}- periBout(run).Tstart(b), periBout(run).fluor{b}(:,fiber(f).ROI(roi) ) ); hold on;
                                %plot( Twindow, predict(fiber(f).(waveType).traceFit{k,roi}, Twindow), 'r', 'lineStyle','--' );
                                xlim([-Inf,Inf]); ylim([-Inf,Inf]); % ylim([-0.02, 1.02]);
                                xlabel('Peri-run time (s)'); ylabel('Normalized Fluor');
                                title(sprintf('[roi, k] = [%i, %i]: effect = %2.2f, R^2 = %2.2f: onset = %2.1f, tauOn = %2.2f, offset = %2.1f, tauOff = %2.2f.  good = %i', ...
                                    roi, k, fiber(f).(waveType).effect(k,roi), fitResult.R2, fitResult.Tonset, fitResult.tauOn, fitResult.Toffset, fitResult.tauOff, fiber(f).(waveType).goodFit(k,roi) ) );
                                
                                subtightplot(3,1,3,opt{:}); cla;
                                plot( periBout(run).T{b} - periBout(run).Tstart(b), periBout(run).velocity{b} );
                                xlim([-Inf,Inf]);
                                xlabel('Peri-run time (s)'); ylabel('Wheel Velocity (cm/s)');
                                pause(0.5);
                            end
                            %}
                        else
                            fprintf('\n[f, roi, run, bout] = [%i, %i, %i, %i] not usable:  %i NaN values', f, roi, run, b, Nnan);
                            fiber(f).(waveType).R2(k,roi) = NaN;
                            fiber(f).(waveType).tauOn(k,roi) = NaN;
                            fiber(f).(waveType).Tonset(k,roi) = NaN;
                        end
                    end
                else
                    fprintf('\nfiber %i, bout %i below minimum length', f, b)
                end
            end
        end
        
        % Fit pairwise separation vs delay to a linear function
        tempROI = find(~all(isnan(fiber(f).(waveType).Tonset),1));
        fiber(f).(waveType).sepDelay = nan(0, 2);
        if numel(tempROI) > 1
            tempPair = nchoosek(tempROI, 2);
            for p = 1:size(tempPair,1)
                pairSep = fiber(f).footSep(tempPair(p,1), tempPair(p,2));
                pairDelays = diff(fiber(f).(waveType).Tonset(:,[tempPair(p,1), tempPair(p,2)]),1,2);
                fiber(f).(waveType).sepDelay = [fiber(f).(waveType).sepDelay; [repmat(pairSep,size(pairDelays,1),1), pairDelays]];
            end
            fiber(f).(waveType).sepDelay(isnan(fiber(f).(waveType).sepDelay(:,2)),:) = [];
        end
        if sum(~isnan(fiber(f).(waveType).sepDelay(:,2))) > 1
            fiber(f).(waveType).sepDelayFit = fitlm(fiber(f).(waveType).sepDelay(:,1), fiber(f).(waveType).sepDelay(:,2), 'RobustOpts',false);
            if fiber(f).(waveType).sepDelayFit.Coefficients.pValue(1) < 0.05
                fiber(f).(waveType).sepDelayConst = fiber(f).(waveType).sepDelayFit.Coefficients.Estimate(1);
            else
                fiber(f).(waveType).sepDelayConst = 0;
            end
            if fiber(f).(waveType).sepDelayFit.Coefficients.pValue(2) < 0.05
                fiber(f).(waveType).sepDelaySlope = fiber(f).(waveType).sepDelayFit.Coefficients.Estimate(2);
            else
                fiber(f).(waveType).sepDelaySlope = 0;
            end
            fiber(f).(waveType).sepDelayPred = polyval([fiber(f).(waveType).sepDelaySlope; fiber(f).(waveType).sepDelayConst], fiber(f).(waveType).sepDelay(:,1));
        else
            fiber(f).(waveType).sepDelayFit = [];
            fiber(f).(waveType).sepDelayConst = NaN;
            fiber(f).(waveType).sepDelaySlope = NaN;
            fiber(f).(waveType).sepDelayPred = [];
        end
        if show && ~isempty(fiber(f).(waveType).sepDelayFit)
            clf;
            plot(fiber(f).(waveType).sepDelay(:,1), fiber(f).(waveType).sepDelay(:,2), '.'); hold on;
            plot( fiber(f).(waveType).sepDelay(:,1), fiber(f).(waveType).sepDelayPred, 'r--' )
            xlabel('Pairwise Separation (um)'); ylabel('Mean Pairwise Lag (s)');
            title(sprintf('Spontaneous Events: f = %i: ', f));
            axis square;
            pause;
        end
    end
    
end
end