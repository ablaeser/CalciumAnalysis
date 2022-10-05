function [waveSpeed, latencyModel, columnGroup, rowGroup, planeGroup ] = GetWaveSpeed(expt, catInfo, periBout, voxInd, varargin) % , waveSpeedPval
% Track how fast a CSD wave washes across a set of specific volumes (ROIs, axons or neuropil)
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
addRequired( IP, 'periBout', @isstruct )
addRequired( IP, 'voxInd', @iscell )
addParameter( IP, 'flatten', true, @islogical ); % [], @isnumeric 
addParameter( IP, 'window', [], @isnumeric) % time window (seconds) in which to detect onset of wave
addParameter( IP, 'method', 'median', @ischar)  % use 'median' method, z threshold crossing ('zCross'), fraction of peak ('fracCross'), or peak of derivative ('deriv')?
addParameter( IP, 'peakFrac', 0.5, @isnumeric )
addParameter( IP, 'zThresh', 3, @isnumeric)
addParameter( IP, 'fitDim', 2, @isnumeric)
addParameter( IP, 'pThresh', 0.01, @isnumeric)
addParameter( IP, 'zoom', 2, @isnumeric) % assume digital zoom = 2
addParameter( IP, 'show', false, @islogical)
addParameter( IP, 'save', '', @ischar)
addParameter( IP, 'saveFig', false, @islogical)
addParameter( IP, 'wait', false, @islogical)
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, catInfo, periBout, voxInd, varargin{:} );  %
flatten = IP.Results.flatten;
onsetWindow = IP.Results.window;
onsetMethod = IP.Results.method; 
pThresh = IP.Results.pThresh;
fitDim = IP.Results.fitDim;
zThresh = IP.Results.zThresh;
peakFrac = IP.Results.peakFrac;
digiZoom = IP.Results.zoom;
show = IP.Results.show;
waitToggle = IP.Results.wait;
saveName = IP.Results.save;
saveFig = IP.Results.saveFig;
overwrite = IP.Results.overwrite;
if flatten, fitDim = 2; end
if ~isempty(saveName),  savePath = sprintf('%s%s_wave_%s.mat', expt.dir, expt.name, saveName);  end
Nvol = numel(voxInd);
latencyModel = cell(periBout.Nbout,Nvol); waveSpeed = nan(periBout.Nbout,Nvol,3,2); %waveSpeedPval = nan(periBout.Nbout,Nvol,3);
gaussFilt = MakeGaussFilt( 4, 0, 2/expt.scanRate, expt.scanRate, false ); %  gaussWidth, gaussMean, gaussSigma, expt(x).scanRate, show 
umPerPixel = (1/0.53)/digiZoom; 
if isempty(saveName) || ~exist(savePath,'file') || overwrite
    for b = 1:periBout.Nbout
        % Get periBout timing and define window (if not already defined)
        Tperi = periBout.T{b} - periBout.Tstart(b);
        XtickZero = find(Tperi == 0);
        XtickPeri = 1:XtickZero-1:numel(Tperi);  
        if isempty( onsetWindow ), onsetWindow = [0, 0.67*periBout.dur(b)]; end
        onsetWindowScans = find(Tperi >= onsetWindow(1) & Tperi <= onsetWindow(2));
        % Get peri-bout volume scans
        periBoutScans = sum(expt.Nscan(1:expt.csd-1)) + periBout.scan{b};
        periBoutVol = double(readSBX(expt.sbx, catInfo, expt.Nplane*(periBoutScans(1)-1)+1, numel(periBoutScans)*expt.Nplane, 1, [])); % double(pipe.io.sbxRead(expt.sbx, expt.Nplane*(periCSDscans(1)-1)+1, numel(periCSDscans)*expt.Nplane, 1, [])); %zeros( )
        periBoutVol = reshape( periBoutVol, expt.Nrow, expt.Ncol, expt.Nplane, periBout.Nscan(b) );
        tic;
        for r = 1:Nvol
            fprintf('\n[b,r] = [%i, %i]  ',b, r);
            if flatten %isempty(zFlat)
                fprintf('Flattening the data...  ');
                tic
                nanScan = nan(expt.Nrow, expt.Ncol, expt.Nplane);
                preFlatVol = nan(expt.Nrow, expt.Ncol, expt.Nplane, periBout.Nscan(b));
                for s = flip(1:periBout.Nscan(b))
                    tempPeriScan = periBoutVol(:,:,:,s);
                    tempPreFlatScan = nanScan;
                    tempPreFlatScan(voxInd{r}) = tempPeriScan(voxInd{r});
                    preFlatVol(:,:,:,s) = tempPreFlatScan;
                end

                % Flatten the volume
                flatVol = squeeze(mean( preFlatVol, 3, 'omitnan'));
                
                %flatVol(isnan(flatVol)) = 0;
                %saveastiff(uint16(flatVol), 'D:\2photon\DL89\171122\testAxon.tif')
                
                firstFrame = flatVol(:,:,1);
                flatInd = find(~isnan(firstFrame(:)));
                [flatRow, flatCol] = ind2sub(size(firstFrame), flatInd);
                Npix = numel(flatInd); %sum(sum(~isnan(flatVol(:,:,1))));
                traceXYZ{b,r} = [flatCol, flatRow, ones(Npix,1)];
                [~,sortInd] = sort( traceXYZ{b,r}(:,1), 'ascend' );
                traceXYZ{b,r} = traceXYZ{b,r}(sortInd,:); 
                % Get traces from flattened volume
                Ntrace(b,r) = Npix;
                rawTrace{b,r} = nan(periBout.Nscan(b), Npix);
                for p = flip(1:Npix)
                    rawTrace{b,r}(:,p) = flatVol(traceXYZ{b,r}(p,2), traceXYZ{b,r}(p,1), :);
                end  
            else
                % Get volume data
                [voxRow, voxCol, voxPlane] = ind2sub([expt.Nrow, expt.Ncol, expt.Nplane], voxInd{r} );
                traceXYZ{b,r} = [voxCol, voxRow, voxPlane]; 
                [~,sortInd] = sort( traceXYZ{b,r}(:,1), 'ascend' );
                traceXYZ{b,r} = traceXYZ{b,r}(sortInd,:); 
                Ntrace(b,r) = size(traceXYZ{b,r},1);
                % Get traces for each voxel
                rawTrace{b,r} = nan(periBout.Nscan(b), Ntrace(b,r));
                for v = flip(1:Ntrace(b,r))
                    rawTrace{b,r}(:,v) = periBoutVol(traceXYZ{b,r}(v,2), traceXYZ{b,r}(v,1), traceXYZ{b,r}(v,3), :);
                end   
            end
            
            % Filter and normalize traces 
            filtTrace{b,r} = filtfilt(gaussFilt, 1, rawTrace{b,r}); %figure; subplot(2,1,1); imagesc(voxrawTrace{b,r}' ); subplot(2,1,2); imagesc(voxfiltTrace{b,r}' ); 
            baseMean{b,r} = mean(filtTrace{b,r}(periBout.preScan{b},:), 1);
            baseStd{b,r} = std(filtTrace{b,r}(periBout.preScan{b},:), 0, 1);
            normTrace{b,r} = (filtTrace{b,r} - baseMean{b,r})./baseStd{b,r}; 
            normDiff{b,r} = [nan(1, Ntrace(b,r)); diff(normTrace{b,r}, 1, 1)];
            rawTraceMean{b,r} = mean( rawTrace{b,r}, 2, 'omitnan' );
            filtTraceMean{b,r} = filtfilt(gaussFilt, 1, rawTraceMean{b,r});
            toc

            % find the wave onset scan for each trace
            fprintf('\nGetting wave onset latencies (%s method)...', onsetMethod);
            traceOnsetScan{b,r} = nan(1,Ntrace(b,r));
            tCheck = find(any( normTrace{b,r}(onsetWindowScans,:) > 3 )); %1:Ntrace(b,r)
            if strcmpi(onsetMethod, 'xcorr')
                % Bin traces into x-y boxes
                binWidth = 10;
                minX = min(traceXYZ{b,r}(:,1)); maxX = max(traceXYZ{b,r}(:,1));
                minY = min(traceXYZ{b,r}(:,2)); maxY = max(traceXYZ{b,r}(:,2));
                xBinLims = minX:binWidth:maxX;
                yBinLims = minY:binWidth:maxY;
                k = 0;
                tic;
                for x = 1:numel(xBinLims)-1
                    for y = 1:numel(yBinLims)-1
                        k = k + 1;
                        tempTraceInd = find( traceXYZ{b,r}(:,1) >= xBinLims(x) & traceXYZ{b,r}(:,1) < xBinLims(x+1) & traceXYZ{b,r}(:,2) >= yBinLims(y) & traceXYZ{b,r}(:,2) < yBinLims(y+1) );
                        binCent(k,:) = mean(traceXYZ{b,r}(tempTraceInd,:));
                        tempCombinedTrace(:,k) = mean( rawTrace{b,r}(:,tempTraceInd), 2, 'omitnan' );
                    end
                    toc
                end
            
            elseif strcmpi(onsetMethod,'fracCross')
                tracePeak = max(normTrace{b,r}(onsetWindowScans,:));
                traceFracPeak = peakFrac*tracePeak;
                if waitToggle, w = waitbar(0, 'Getting onsets'); end
                %figure;
                for t = tCheck
                    try  %#ok<*TRYNC>
                        traceOnsetScan{b,r}(t) = find(normTrace{b,r}(onsetWindowScans,t) > traceFracPeak(t), 1,'first' ); 
                        %{ 
                        plot( normTrace{b,r}(:,t) ); hold on;
                        plot(onsetWindowScans, normTrace{b,r}(onsetWindowScans,t), 'k' );
                        k = traceOnsetScan{b,r}(t) + onsetWindowScans(1) - 1;
                        plot( k, normTrace{b,r}(k,t), 'o');
                        pause; cla;
                        %}
                    end
                    if waitToggle, waitbar(t/Ntrace(b,r), w); end
                end
                if waitToggle, delete(w); end
            elseif strcmpi(onsetMethod,'median')
                if waitToggle, w = waitbar(0, 'Getting onsets'); end
                %figure;
                for t = tCheck 
                    try 
                        baseQuartiles = prctile( normTrace{b,r}(periBout.preScan{b},t), [25,50,75]);
                        medianRuleThresh = baseQuartiles(2) + 2.3*diff(baseQuartiles([1,3]));
                        traceOnsetScan{b,r}(t) = find(normTrace{b,r}(onsetWindowScans,t) >= medianRuleThresh, 1,'first' ); 
                        %{ 
                        plot( normTrace{b,r}(:,t) ); hold on;
                        plot(onsetWindowScans, normTrace{b,r}(onsetWindowScans,t), 'k' );
                        k = traceOnsetScan{b,r}(t) + onsetWindowScans(1) - 1;
                        plot( k, normTrace{b,r}(k,t), 'o');
                        pause; cla;
                        %}
                    end
                    %if waitToggle, waitbar(t/Ntrace(b,r), w); end
                end
                if waitToggle, delete(w); end
            elseif strcmpi(onsetMethod,'xCorr')
                
                %{
                X = rawTraceMean{b,r};
                for t = tCheck
                    Y = rawTrace{b,r}(:,t);
                    [traceXcorr, lags] = xcorr( X, Y ); %xcorr( normTrace{b,r}(:,1), normTrace{b,r}(:,end) );
                    [XcorrPeak, XcorrPeakLag] = max(traceXcorr);
                    subplot(2,1,1); cla; plot( [X,Y] );
                    subplot(2,1,2); cla; plot(lags, traceXcorr ); hold on;
                    plot( lags(XcorrPeakLag), XcorrPeak, 'o');
                    pause(0.1);
                end
                
                % Bin data by groups of columns (percentile)
                groupPct = 10;
                Ngroup = ceil(100/groupPct); 
                for g = flip(1:Ngroup)
                    stimGroup.lim(g,:) = [prctile(stim,(g-1)*groupPct), prctile(stim,g*groupPct)];
                    stimGroup.ind{g} = find( stim >= prctile(stim,groupPct*(g-1)) & stim < prctile(stim,groupPct*g) ); 
                end
                %}
            elseif strcmpi(onsetMethod,'zCross')
                w = waitbar(0, 'Getting onsets');
                for t = tCheck 
                    try
                        traceOnsetScan{b,r}(t) = find(normTrace{b,r}(onsetWindowScans,t) > zThresh, 1,'first' ); 
                        waitbar(t/Ntrace(b,r), w);
                    end
                end
                delete(w);
            elseif strcmpi(onsetMethod,'deriv')
                [~, traceOnsetScan{b,r}(tCheck)] = max( normDiff{b,r}(onsetWindowScans,tCheck), [], 1 ); 
            end
            if all(isnan(traceOnsetScan{b,r})), error('No good onsets found'); end
            traceOnsetScan{b,r} = traceOnsetScan{b,r} + onsetWindowScans(1) - 1; % correct for window offset
            goodTrace = find( ~isnan(traceOnsetScan{b,r}) & baseMean{b,r} > 0 ); % Suppress voxels that are all zeros, or where no onset was detected
            NgoodTrace(b,r) = numel(goodTrace);
            onsetLatency{b,r} = nan(1,Ntrace(b,r));
            onsetLatency{b,r}(goodTrace) = Tperi(traceOnsetScan{b,r}(goodTrace));
            
            % Form a 3D latency volume and fit a linear model
            if flatten
                latencyVolume{b,r} = nan(expt.Nrow, expt.Ncol, 1);
            else
                latencyVolume{b,r} = nan(expt.Nrow, expt.Ncol, expt.Nplane);
            end
            for t = goodTrace%1:Ntrace(b,r)
                latencyVolume{b,r}(traceXYZ{b,r}(t,2), traceXYZ{b,r}(t,1), traceXYZ{b,r}(t,3)) = onsetLatency{b,r}(t);
            end
            %imagesc( latencyVolume)
            if fitDim == 2 || fitDim == 3
                fprintf('\nFitting the %iD model...', fitDim);
                latencyModel{b,r} = fitlm( traceXYZ{b,r}(:,1:fitDim), onsetLatency{b,r}', 'RobustOpts','on' );
                latencyCoeff = latencyModel{b,r}.Coefficients; % seconds per pixel or plane
                % Convert coefficients to speeds
                conversionFactor = [umPerPixel, umPerPixel, 1];
                for d = 1:fitDim
                    waveSpeed(b,r,d,2) = latencyCoeff.pValue(d+1);
                    if waveSpeed(b,r,d,2) < pThresh
                        waveSpeed(b,r,d,1) = conversionFactor(d)*(1/latencyCoeff.Estimate(d+1)); 
                    end
                end
            else
                error('fitDim must be 2 or 3');
            end

            % Group latencies by column, row and bin
            columnGroup(b,r) = struct('N',0, 'value',[], 'latency',[], 'tickInt',30, 'ticks',[], 'ticksLabels',[]);
            columns = unique( traceXYZ{b,r}(:,1)); % , ~, columnInd    , 'stable'
            columnGroup(b,r).N = numel(columns);
            columnGroup(b,r).value = columns';
            for p = 1:columnGroup(b,r).N
                columnGroup(b,r).latency{p} = onsetLatency{b,r}( traceXYZ{b,r}(:,1) == columns(p) );
            end
            if columnGroup(b,r).N < columnGroup(b,r).tickInt, columnGroup(b,r).tickInt = 5; end
            columnGroup(b,r).ticks = 1:columnGroup(b,r).tickInt:columnGroup(b,r).N; 
            columnGroup(b,r).ticksLabels = columnGroup(b,r).value(columnGroup(b,r).ticks);

            rowGroup(b,r) = struct('N',0, 'value',[], 'latency',[], 'tickInt',30, 'ticks',[], 'ticksLabels',[]);
            rows = unique( traceXYZ{b,r}(:,2)); % , ~, rowInd    , 'stable'
            rowGroup(b,r).N = numel(rows);
            rowGroup(b,r).value = rows';
            for p = 1:rowGroup(b,r).N
                rowGroup(b,r).latency{p} = onsetLatency{b,r}( traceXYZ{b,r}(:,2) == rows(p) );
            end
            if rowGroup(b,r).N < rowGroup(b,r).tickInt, rowGroup(b,r).tickInt = 5; end
            rowGroup(b,r).ticks = 1:rowGroup(b,r).tickInt:rowGroup(b,r).N; 
            rowGroup(b,r).ticksLabels = rowGroup(b,r).value(rowGroup(b,r).ticks);
            
            planeGroup(b,r) = struct('N',0, 'value',[], 'latency',[], 'tickInt',3, 'ticks',[], 'ticksLabels',[]);
            if fitDim == 3
                planes = unique( traceXYZ{b,r}(:,3)); % , ~, planeInd    , 'stable'
                planeGroup(b,r).N = numel(planes);
                planeGroup(b,r).value = planes';
                for p = 1:planeGroup(b,r).N
                    planeGroup(b,r).latency{p} = onsetLatency{b,r}( traceXYZ{b,r}(:,3) == planes(p) );
                end
                planeGroup(b,r).ticks = 1:planeGroup(b,r).tickInt:planeGroup(b,r).N; 
                planeGroup(b,r).ticksLabels = planeGroup(b,r).value(planeGroup(b,r).ticks);
            end
            toc
        end
    end
    clearvars flatVol flatRow flatCol flatInd preFlatVol firstFrame nanScan voxCol voxRow voxPlane tempPreFlatScan tempPeriScan periBoutVol 
    if ~isempty(saveName)
        fprintf('\nSaving %s', savePath); % , 'normTrace','normDiff','Tperi','onsetWindow','onsetMethod','
        save(savePath, '-v7.3');   
    end
else
    fprintf('\nLoading %s', savePath);
    load(savePath, 'latencyModel', 'waveSpeed', 'columnGroup', 'rowGroup', 'planeGroup'); 
end
toc;

% Plot the process (optional)
if show
    opt = {[0.1,0.07], [0.1,0.04], [0.05,0.05]};  % {[vert, horz], [bottom, top], [left, right] }
    close all; 
    WaveSpeedFig = figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
    for b = 1:periBout.Nbout
        for r = 1:Nvol
            medLatRange = [0, Inf]; % prctile(latencyVolume{b,r}(:), [3,97] prctile(latencyVolume{b,r}(:), 97)
            subtightplot(3,3,[1,4,7],opt{:}); cla; % sp(1) =
            if strcmpi(onsetMethod,'median')    
                imagesc(normTrace{b,r}'); hold on;
                CB(1) = colorbar; CB(1).Label.String = 'z-score'; CB(1).Label.FontWeight = 'bold';  
                caxis([-1,5]);
                title( sprintf('Onsets detected for %i of %i traces (median method)', NgoodTrace(b,r), Ntrace(b,r) ), 'FontSize',10 );
            elseif strcmpi(onsetMethod,'zCross')
                imagesc(normTrace{b,r}'); hold on;
                CB(1) = colorbar; CB(1).Label.String = 'z-score'; CB(1).Label.FontWeight = 'bold';  
                caxis([-1,5]);
                title( sprintf('Onsets detected for %i of %i traces (z>%2.1f)', NgoodTrace(b,r), Ntrace(b,r), zThresh ), 'FontSize',10 );
            elseif strcmpi(onsetMethod,'fracCross')
                imagesc(normTrace{b,r}'); hold on;
                CB(1) = colorbar; CB(1).Label.String = 'z-score'; CB(1).Label.FontWeight = 'bold';  
                caxis([-1,5]);
                title( sprintf('Onsets detected for %i of %i traces (cross %2.f pct of peak)', NgoodTrace(b,r), Ntrace(b,r), 100*peakFrac ), 'FontSize',10 );
            elseif strcmpi(onsetMethod,'deriv')
                imagesc(normalize(normDiff{b,r},1)'); hold on;
                CB(1) = colorbar; CB(1).Label.String = 'dF/dt'; CB(1).Label.FontWeight = 'bold'; 
                title( sprintf('Onsets detected for %i of %i traces (max derivative)', NgoodTrace(b,r), Ntrace(b,r) ), 'FontSize',10 ); %title( 'Normalized Voxel Derivatives' );
            end
            if Ntrace(b,r) < 5000, plot( traceOnsetScan{b,r}, 1:Ntrace(b,r), 'k.', 'markersize',1); end
            line( (onsetWindowScans(1)-0.5)*[1,1], [1, Ntrace(b,r)], 'color','r', 'LineStyle','--');
            line( (onsetWindowScans(end)+0.5)*[1,1], [1, Ntrace(b,r)], 'color','r', 'LineStyle','--');
            set(gca, 'Xtick',XtickPeri, 'XtickLabel',sprintfc('%2.1f', Tperi(XtickPeri)), 'Ytick',[] );
            ylabel('< Right  Trace  Left >'); xlabel('Peri-Wave Time (s)');
            if flatten
                subtightplot(3,3,[2,3,5,6],opt{:}); cla;
                imagesc( latencyVolume{b,r} );
                axis image;
                CB(2) = colorbar('Location','EastOutside'); CB(2).Label.String = 'Onset Latency (s)'; CB(2).Label.FontWeight = 'bold'; 
                caxis(medLatRange) %caxis([min(medLatencyFoot(:)), max(medLatencyFoot(:))]) % onsetWindow prctile(latencyVolume{b,r}(:), [1,99]
                set(gca,'Xtick',[], 'Ytick',[]);
                xlabel('AP'); ylabel('ML');
                title( sprintf('%s: Bout %i, Volume %i', expt.name, b, r), 'Interpreter','none' ); 
            else
                subtightplot(3,3,2,opt{:}); cla;
                imagesc( median(latencyVolume{b,r}, 3, 'omitnan') ); axis image;
                CB(2) = colorbar; CB(2).Label.String = 'Median Onset Latency (s)'; CB(2).Label.FontWeight = 'bold'; 
                caxis(medLatRange) %caxis([min(medLatencyFoot(:)), max(medLatencyFoot(:))]) % onsetWindow prctile(latencyVolume{b,r}(:), [1,99]
                set(gca,'Xtick',[], 'Ytick',[]);
                xlabel('AP'); ylabel('ML');
                title( sprintf('%s: Volume %i', expt.name, r), 'Interpreter','none' ); 

                subtightplot(3,3,5,opt{:}); cla;
                imagesc( squeeze(median(latencyVolume{b,r}, 2, 'omitnan'))' ); axis image;
                set(gca,'Xtick',[], 'Ytick',[]);
                caxis(medLatRange) 
                xlabel('ML'); ylabel('Z');

                subtightplot(3,3,8,opt{:}); cla;
                imagesc( squeeze(median(latencyVolume{b,r}, 1, 'omitnan'))' ); axis image;
                set(gca,'Xtick',[], 'Ytick',[]);
                caxis(medLatRange) 
                xlabel('AP'); ylabel('Z');
            end

            if flatten || fitDim == 2
                subtightplot(3,3,8,opt{:}); cla;
                boxplot( cell2padmat(columnGroup(b,r).latency) );
                ylim(medLatRange);
                set(gca, 'Xtick',columnGroup(b,r).ticks, 'XtickLabel', columnGroup(b,r).ticksLabels, 'FontSize',8)
                xlabel('AP Position (column)'); ylabel('Latency (s)'); xtickangle(30);
                title(sprintf('Estimated speed component = %2.2f (um/s)  (p=%0.4f)', waveSpeed(b,r,1,1), waveSpeed(b,r,1,2)) );

                subtightplot(3,3,9,opt{:}); cla;
                boxplot( cell2padmat(rowGroup(b,r).latency) );
                ylim(medLatRange);
                set(gca, 'Xtick',rowGroup(b,r).ticks, 'XtickLabel', rowGroup(b,r).ticksLabels, 'FontSize',8)
                xlabel('ML Position (row)'); ylabel('Latency (s)'); xtickangle(30);
                title(sprintf('Estimated speed component = %2.2f  (p=%0.4f)', waveSpeed(b,r,2,1), waveSpeed(b,r,2,2)) );
            else
                subtightplot(3,3,3,opt{:}); cla;
                boxplot( cell2padmat(columnGroup(b,r).latency) );
                ylim(medLatRange);
                set(gca, 'Xtick',columnGroup(b,r).ticks, 'XtickLabel', columnGroup(b,r).ticksLabels)
                xlabel('AP Position (column)'); ylabel('Latency (s)'); 
                title(sprintf('Estimated speed component = %2.2f (um/s)  (p=%0.4f)', waveSpeed(b,r,1,1), waveSpeed(b,r,1,2)) );

                subtightplot(3,3,6,opt{:}); cla;
                boxplot( cell2padmat(rowGroup(b,r).latency) );
                ylim(medLatRange);
                set(gca, 'Xtick',rowGroup(b,r).ticks, 'XtickLabel', rowGroup(b,r).ticksLabels, 'FontSize',8)
                xlabel('ML Position (row)'); ylabel('Latency (s)');
                xtickangle(30);
                title(sprintf('Estimated speed component = %2.2f  (p=%0.4f)', waveSpeed(b,r,2,1), waveSpeed(b,r,2,2)) );

                subtightplot(3,3,9,opt{:}); cla;
                boxplot( cell2padmat(planeGroup(b,r).latency) );
                ylim(medLatRange);
                set(gca, 'Xtick',planeGroup(b,r).ticks, 'XtickLabel', planeGroup(b,r).ticksLabels, 'FontSize',8)

                xlabel('Z Position (plane)'); ylabel('Latency (s)');
                title(sprintf('Estimated speed component = %2.2f  (p=%0.4f)', waveSpeed(b,r,3,1), waveSpeed(b,r,3,2)) );
                %if periBout.Nbout > 1 || Nvol > 1,  pause;   end  %(1);
            end

            if saveFig
                figPath = sprintf('%s%s_WaveSpeed_%s_bout%i_vol%i.tif', expt.dir, expt.name, saveName, b, r);
                fprintf('\nSaving %s', figPath);
                print(WaveSpeedFig, figPath, '-dtiff')
            else
                impixelinfo;
                pause;
            end
        end
    end
end
end