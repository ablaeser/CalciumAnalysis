function csdWave = GetNeuropilWaveSpeed(expt, catInfo, periBout, fluor, voxInd, varargin) % , waveSpeedPval
% Track how fast a CSD wave washes across a set of specific volumes (ROIs, axons or neuropil)
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
addRequired( IP, 'periBout', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'voxInd', @iscell )
addParameter( IP, 'flatten', true, @islogical ); % [], @isnumeric 
addParameter( IP, 'bin', 50, @isnumeric) % width of binning (pixels)
addParameter( IP, 'window', [], @isnumeric) % time window (seconds) in which to detect onset of wave
addParameter( IP, 'zThresh', 3, @isnumeric)
addParameter( IP, 'fitDim', 2, @isnumeric)
addParameter( IP, 'pThresh', 0.01, @isnumeric)
addParameter( IP, 'minXC', 0.5, @isnumeric )
addParameter( IP, 'minSep', 100/expt.scanRate, @isnumeric ) % 150
addParameter( IP, 'maxDelay', 10, @isnumeric )
addParameter( IP, 'show', false, @islogical)
addParameter( IP, 'save', '', @ischar)
addParameter( IP, 'saveFig', false, @islogical)
addParameter( IP, 'wait', false, @islogical)
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, catInfo, periBout, fluor, voxInd, varargin{:} );  %
flatten = IP.Results.flatten;
onsetWindow = IP.Results.window;
binWidth = IP.Results.bin;
%onsetMethod = IP.Results.method; 
%pThresh = IP.Results.pThresh;
%fitDim = IP.Results.fitDim;
zThresh = IP.Results.zThresh;
minXC = IP.Results.minXC;
minSep = IP.Results.minSep;
maxDelay = IP.Results.maxDelay;
%peakFrac = IP.Results.peakFrac;
show = IP.Results.show;
%waitToggle = IP.Results.wait;
saveName = IP.Results.save;
%saveFig = IP.Results.saveFig;
overwrite = IP.Results.overwrite;

if ~isempty(saveName),  savePath = sprintf('%s%s_wave_%s.mat', expt.dir, expt.name, saveName);  end
%Nvol = numel(voxInd);
gaussFilt = MakeGaussFilt( 4, 0, 1/expt.scanRate, expt.scanRate, false ); %  gaussWidth, gaussMean, gaussSigma, expt(x).scanRate, show 

%close all; clearvars sp;  
figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
if isempty(saveName) || ~exist(savePath,'file') || overwrite
    for b = 1:periBout.Nbout
        % Get peri-bout volume scans
        periBoutScans = sum(expt.Nscan(1:expt.csd-1)) + periBout.scan{b};
        periBoutVol = double(readSBX(expt.sbx, catInfo, expt.Nplane*(periBoutScans(1)-1)+1, numel(periBoutScans)*expt.Nplane, 1, [])); % double(pipe.io.sbxRead(expt.sbx, expt.Nplane*(periCSDscans(1)-1)+1, numel(periCSDscans)*expt.Nplane, 1, [])); %zeros( )
        periBoutVol = reshape( periBoutVol, expt.Nrow, expt.Ncol, expt.Nplane, periBout.Nscan(b) );

        % Flatten the volume
        fprintf('Flattening the data...  ');
        tic
        nanScan = nan(expt.Nrow, expt.Ncol, expt.Nplane);
        preFlatVol = nan(expt.Nrow, expt.Ncol, expt.Nplane, periBout.Nscan(b));
        for s = flip(1:periBout.Nscan(b))
            tempPeriScan = periBoutVol(:,:,:,s);
            tempPreFlatScan = nanScan;
            tempPreFlatScan(voxInd{1}) = tempPeriScan(voxInd{1});
            preFlatVol(:,:,:,s) = tempPreFlatScan;
        end
        flatVol = squeeze(mean( preFlatVol, 3, 'omitnan'));
        %flatVol(isnan(flatVol)) = 0;
        %saveastiff(uint16(flatVol), 'D:\2photon\DL89\171122\testAxon.tif')
        firstFrame = flatVol(:,:,1);
        flatInd = find(~isnan(firstFrame(:)));
        [flatRow, flatCol] = ind2sub(size(firstFrame), flatInd);
        Npix = numel(flatInd); %sum(sum(~isnan(flatVol(:,:,1))));
        traceXYZ{b} = [flatCol, flatRow, ones(Npix,1)];
        [~,sortInd] = sort( traceXYZ{b}(:,1), 'ascend' );
        [~,ySortInd] = sort( traceXYZ{b}(:,2), 'ascend' );
        traceXYZ{b} = traceXYZ{b}(sortInd,:); 
        % Get traces from flattened volume
        rawTrace{b} = nan(periBout.Nscan(b), Npix);
        for p = flip(1:Npix)
            rawTrace{b}(:,p) = flatVol(traceXYZ{b}(p,2), traceXYZ{b}(p,1), :);
        end  
        toc

        % Get periBout timing and define window (if not already defined)
        onsetWindow = [-10,30];
        Tperi = periBout.T{b} - periBout.Tstart(b);
        if isempty( onsetWindow ), onsetWindow = [-10, 0.67*periBout.dur(b)]; end
        windowScans = find(Tperi >= onsetWindow(1) & Tperi <= onsetWindow(2));
        Twindow = Tperi(windowScans);
        % Bin traces into x-y boxes
        tic
        fprintf('\nBinning pixels (binWidth = %i)', binWidth)
        minX = min(traceXYZ{b}(:,1)); maxX = max(traceXYZ{b}(:,1));
        minY = min(traceXYZ{b}(:,2)); maxY = max(traceXYZ{b}(:,2));
        xBinLims = minX:binWidth:maxX;
        yBinLims = minY:binWidth:maxY;
        NbinX = numel(xBinLims)-1; NbinY = numel(yBinLims)-1;
        Nbin = NbinX*NbinY; %size(binTrace,2);
        binTrace = nan(NbinX, NbinY, periBout.Nscan(b)); binCent = nan(NbinX, NbinY, 3);
        for x = 1:NbinX
            for y = 1:NbinY
                boxTraceInd = find( traceXYZ{b}(:,1) >= xBinLims(x) & traceXYZ{b}(:,1) < xBinLims(x+1) & traceXYZ{b}(:,2) >= yBinLims(y) & traceXYZ{b}(:,2) < yBinLims(y+1) );
                binTrace(x,y,:) = filtfilt(gaussFilt, 1, mean(rawTrace{b}(:,boxTraceInd), 2, 'omitnan') ); % filtTrace
                binCent(x,y,:) = expt.umPerPixel*mean(traceXYZ{b}(boxTraceInd,:));
            end
        end
        %saveastiff( uint16(permute(binTrace, [2,1,3])), 'D:\2photon\DL115\180707\binTrace.tif' )
        %preScan = find(Tperi < 0);
        %binPreMean = mean( binTrace(:,:,preScan), 3 );
        %binPreStd = std( binTrace(:,:,preScan), 0, 3 );
        %binTrace = (binTrace - binPreMean)./binPreStd;  %min( binTrace(:,:,Tperi<0),[], 3 )
        %binTrace = normalize( binTrace, 3); % (:,:,onsetWindowScans)
        binTrace = (binTrace - min(binTrace,[],3)); %./max(binTrace,[],3)
        binTrace = binTrace./max(binTrace,[],3);
        binTrace = binTrace(:,:,windowScans);
        
        %binTrace = binTrace - min( binTrace(:,:,Tperi<0),[], 3 );
        %binTrace = binTrace./max(binTrace(:,:,Tperi>0),[], 3 );
        %plot( squeeze(binTrace(1,1,:)) );
        %plot( squeeze(binTrace(NbinX,1,:)) );
        toc
            
        %{ 
        % Estimate delay between each bin and the overall CSD signal, using xcorr
        %csdSignal = fluor.z.vol(onsetWindowScans,:);
        k = 0; csdOnsetSignal = nan(numel(onsetWindowScans), 1);
        for s = onsetWindowScans'
            k = k+1;
            tempCSDscan = flatVol( min(flatRow):max(flatRow), min(flatCol):max(flatCol), s );
            csdOnsetSignal(k) = mean( tempCSDscan(:), 'omitnan' );
        end
        csdOnsetSignal = normalize(csdOnsetSignal,1);
        figure;
        binCSDlag = nan(NbinX, NbinY);
        for x = 1:NbinX
            for y = 1:NbinY
                [tempPairXC, tempPairLagScans] = xcorr( squeeze(binTrace(x,y,:)), squeeze(binTrace(1,1,:)), round(maxDelay/2*expt.scanRate), 'normalized' ); %  csdOnsetSignal
                [binXCpeak, binXClag] = max( tempPairXC );
                if binXCpeak > minXC 
                    binCSDlag(x,y) = tempPairLagScans(binXClag)/expt.scanRate;
                end
                %{
                subplot(1,3,1); plot(  [squeeze(binTrace(x,y,:)), squeeze(binTrace(1,1,:))] );  % csdOnsetSignal
                axis square;
                subplot(1,3,2); plot(tempPairLagScans, tempPairXC ); axis square;
                pause%(0.1);
                %}
            end
        end
        imagesc( binCSDlag' ); axis image; impixelinfo; % subplot(1,3,3); 
        %}
            
        % Estimate delay between pairs of well-separated, highly cross-correlated bins
        Ndelay = round(maxDelay*expt.scanRate);
        delayVec = (-Ndelay:Ndelay)/expt.scanRate;
        dataMat = nan(Nbin^2, 5);
        k = 0;
        tic;
        %figure;
        for x = 1:NbinX-1
            for y = 1:NbinY-1
                for X = x:NbinX %1:NbinX
                    for Y = y:NbinY %1:NbinY
                        tempDX = binCent(X,Y,1) - binCent(x,y,1);
                        tempDY = binCent(X,Y,2) - binCent(x,y,2);
                        tempSep = sqrt(tempDX^2 + tempDY^2 );
                        if tempSep > minSep
                            % x= 1; y = 1; X = NbinX; Y = 1; %NbinY;
                            trace1 = squeeze(binTrace(X,Y,:));
                            trace2 = squeeze(binTrace(x,y,:));
                            tempPairXC = xcorr( trace1,  trace2, Ndelay, 'normalized' );
                            [tempXCpeak, tempPeakDelayInd] = max( tempPairXC );
                            if tempXCpeak > minXC
                                % use parabola interpolation to improve estimate of peak location (https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html)
                                tempPairXCtop3 = tempPairXC(tempPeakDelayInd-1:tempPeakDelayInd+1); % alpha, beta, gamma
                                p = 0.5*(tempPairXCtop3(1)-tempPairXCtop3(3))/(tempPairXCtop3(1)-2*tempPairXCtop3(2)+tempPairXCtop3(3));
                                estPeakDelayInd = tempPeakDelayInd + p;
                                estPeakDelay = interp1(1:2*Ndelay+1, delayVec, estPeakDelayInd );
                                estPeakXC = interp1(1:2*Ndelay+1, tempPairXC, estPeakDelayInd );
                                %{
                                subplot(1,2,1); plot([trace1, trace2] ); title( sprintf('[x,y] = [%i, %i] (red),  [X,Y] = [%i, %i] (blue)', x, y, X, Y) );
                                subplot(1,2,2); plot( delayVec, tempPairXC ); hold on;
                                plot(delayVec(tempPeakDelayInd), tempPairXC(tempPeakDelayInd), 'o' ); 
                                plot( estPeakDelay, estPeakXC, 'kx' );
                                title('Peak (circle) and estimated peak (x)');
                                xlabel('Delay (s, red leads blue)'); ylabel('Cross-Correlation');
                                hold off;
                                pause
                                %}
                                k = k + 1;
                                dataMat(k,:) = [tempDX, tempDY, tempSep, estPeakDelay, estPeakXC]; % , delayVec(tempPeakDelayInd), tempXCpeak
                            end
                        end
                    end
                end
            end
        end
        dataMat(k+1:end,:) = [];
        toc;
        
        dataMdl = fitlm(dataMat(:,1:2), dataMat(:,4), 'RobustOpts',true); % , 'y ~ x1 + x2 - 1'
        csdWave.velocity = 1./dataMdl.Coefficients.Estimate(2:3);
        [csdWave.angle, csdWave.speed] = cart2pol(csdWave.velocity(1), -csdWave.velocity(2)); % y axis is flipped
        
        binTraceXsort = zeros( size(binTrace,3), Nbin);  binTraceYsort = zeros( size(binTrace,3), Nbin);
        k = 0;
        for x = 1:NbinX
            for y = 1:NbinY
                k = k+1;
                binTraceXsort(:,k) = squeeze(binTrace(x,y,:));
            end
        end
        k = 0;
        for y = 1:NbinY
            for x = 1:NbinX
                k = k+1;
                binTraceYsort(:,k) = squeeze(binTrace(x,y,:));
            end
        end

        opt = {[0.10,0.09], [0.08,0.06], [0.08,0.04]};  % {[vert, horz], [bottom, top], [left, right] }
        close all; clearvars sp;  
        XtickZero = find(Twindow == 0);
        XtickWindow = 1:XtickZero-1:numel(Twindow); 
      
        figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
        subplot(2,3,1); 
        imagesc( binTraceXsort' ); %axis image;
        cb = colorbar; cb.Label.String = 'z-score';
        xlabel('Peri-CSD Time (s)');  ylabel('Tile');
        set(gca, 'Ytick',[1,Nbin], 'YtickLabel',{'Leftmost','Rightmost'}, 'Xtick',XtickWindow, 'XtickLabel',sprintfc('%2.1f', Twindow(XtickWindow)), 'TickDir','out' ); % 
        
        subplot(2,3,2); 
        plot(dataMat(:,4), dataMat(:,1), '.' ); hold on;
        
        axis square;
        xlim([-5,5]); ylim([0,800])
        xlabel('dT (s)'); ylabel('dX (um)'); 
        title(sprintf('X axis: Estimated velocity = %2.1f (p = %2.4f)', csdWave.velocity(1), dataMdl.Coefficients.pValue(2))); 
        
        subplot(2,3,4); 
        imagesc( binTraceYsort' ); %axis square;%axis image;
        cb = colorbar; cb.Label.String = 'z-score';
        xlabel('Peri-CSD Time (s)');  ylabel('Tile');
        set(gca, 'Ytick',[1,Nbin], 'YtickLabel',{'Topmost','Bottommost'}, 'Xtick',XtickWindow, 'XtickLabel',sprintfc('%2.1f', Twindow(XtickWindow)), 'TickDir','out' ); 
        xlim([-Inf,Inf]);
        
        subplot(2,3,5); 
        plot(dataMat(:,4), dataMat(:,2), '.' ); 
        xlim([-5,5]); ylim([0,800])
        xlabel('dT (s)'); ylabel('dY (um)'); axis square;
        title(sprintf('Y axis: Estimated velocity = %2.1f (p = %2.4f)', csdWave.velocity(2), dataMdl.Coefficients.pValue(3)));
        
        subplot(2,3,[3,6]);
        polarplot( csdWave.angle*[1,1], csdWave.speed*[0,1], 'LineWidth',3, 'Color','r' ); hold off;

        
        

        
        %{
        fprintf('\nCalculating pairwise cross-correlations')
        tempPairInd = find( triu(binSep,1) > minSep ); %  & XCpeak > minXC
        clearvars tempPair
        [tempPair(:,1), tempPair(:,2)] = ind2sub( [Nbin, Nbin], tempPairInd ); 
        Npair = size(tempPair,1);
        % Find peak xcorr for first pair of bins
        XCpeak = nan(Nbin, Nbin); XCpeakLag = nan(Nbin, Nbin);    
        [tempPairXC, tempPairLagScans] = xcorr( binTrace(:,tempPair(1,1)),  binTrace(:,tempPair(1,2)), round(maxDelay*expt.scanRate), 'normalized' ); 
        pairLags = (1/expt.scanRate)*((1:numel(tempPairLagScans))  - median(1:numel(tempPairLagScans)));
        [XCpeak(tempPair(1,1), tempPair(1,2)), pairXClagScan] = max( tempPairXC );
        XCpeakLag(tempPair(1,1), tempPair(1,2)) = pairLags(pairXClagScan);
        w = waitbar(1/Npair, 'Calculating pairwise cross-correlations');
        for p = 2:Npair
            [tempPairXC, tempPairLagScans] = xcorr( binTrace(:,tempPair(p,1)),  binTrace(:,tempPair(p,2)), round(maxDelay*expt.scanRate), 'normalized' ); %xcorr( binTrace(:,[tempPair(p,1),tempPair(p,2)]), round(maxDelay*expt.scanRate), 'normalized' );
            [XCpeak(tempPair(p,1), tempPair(p,2)), pairXClagScan] = max( tempPairXC ); %XCpeak(p), XCpeakLagScan(p)
            XCpeakLag(tempPair(p,1), tempPair(p,2)) = pairLags(pairXClagScan);
            waitbar(p/Npair, w)
        end
        delete(w);
        toc

        % Estimate CSD conduction velocity
        tempPairInd = find(XCpeak > minXC); %find( triu(binSep,1) > minSep & XCpeak > minXC ); %use only well-separated, highly correlated pairs of ROI
        if ~isempty(tempPairInd)
            clearvars tempPair;
            [tempPair(:,1), tempPair(:,2)] = ind2sub( [Nbin, Nbin], tempPairInd ); 
            csdWave.pairDisp = binCent(tempPair(:,1),1:2) - binCent(tempPair(:,2),1:2); % displacement vectors between pairs
            %csdWave.pairDisp(abs(csdWave.pairDisp) < minSep) = NaN; % suppress short dimensions
            csdWave.pairLag = XCpeakLag(tempPairInd);
            csdWave.pairVelocity = csdWave.pairDisp./repmat(csdWave.pairLag,1,2);
            csdWave.velocity = median(csdWave.pairVelocity, 1, 'omitnan');
            [csdWave.angle, csdWave.speed] = cart2pol(csdWave.velocity(1), -csdWave.velocity(2)); % y axis is flipped
            %csdWave.angle = rad2deg(csdWave.angle);
            %tempMdl = fitlm( csdWave.pairDisp, csdWave.pairLag, 'RobustOpts','on' );
            %tempMdlCoeff = 1./tempMdl.Coefficients.Estimate;

        end
        %}

        if show
            XtickZero = find(Tperi == 0);
            XtickWindow = 1:XtickZero-1:numel(Tperi); 
            opt = {[0.10,0.09], [0.08,0.06], [0.08,0.04]};  % {[vert, horz], [bottom, top], [left, right] }
            % PLOT CSD PROCESS
            subtightplot(2,3,1, opt{:}); 
            imagesc( binTrace'); axis square;
            cb = colorbar; cb.Label.String = 'z-score';
            xlabel('Peri-CSD Time (s)');  ylabel('Tile');
            set(gca, 'Ytick',[1,Nbin], 'YtickLabel',{'Leftmost','Rightmost'}, 'Xtick',XtickWindow, 'XtickLabel',sprintfc('%2.1f', Tperi(XtickWindow)) ); % 
            xlim([-Inf,Inf]);
            title( sprintf('Estimated conduction velocity = [%2.1f,  %2.1f] um/s  (AP, ML)',csdWave.velocity ) ); % ,  %2.1f

            subtightplot(2,3,2, opt{:});
            imagesc( binTrace(:,yBinSort)' );  % 
            axis square;
            xlabel('Peri-CSD Time (s)');  ylabel('Tile');
            set(gca, 'Ytick',[1,Nbin], 'YtickLabel',{'Topmost','Bottommost'}, 'Xtick',XtickWindow, 'XtickLabel',sprintfc('%2.1f', Tperi(XtickWindow)) ); %
            xlim([-Inf,Inf]);

            subtightplot(2,3,3, opt{:});
            % {
            showPairVelocity = csdWave.pairVelocity; showPairVelocity(isnan(showPairVelocity))= 0;
            [pairTheta, pairMag] = cart2pol( showPairVelocity(:,1), showPairVelocity(:,2) );
            for p = 1:length(pairTheta)
                 polarplot( pairTheta(p)*[1,1], pairMag(p)*[0,1], 'LineWidth',1, 'Color',[0,0,0,0.03] ); 
                 hold on; 
                 %pause%(0.01);
            end
            showVelocity = csdWave.velocity; showVelocity(isnan(showVelocity))= 0;
            %}
            [medTheta, medMag] = cart2pol( csdWave.velocity(:,1), csdWave.velocity(:,2) );
            polarplot( medTheta*[1,1], medMag*[0,1], 'LineWidth',3, 'Color','r' ); hold off;

            subtightplot(2,3,4, opt{:});  
            imagesc( binSep ); axis square;
            colorbar;
            xlabel('Tile'); ylabel('Tile'); 
            title( sprintf('Pairwise Distance (min = %2.2f um)', minSep )); %title('Separation'); 
            set(gca,'Xtick',[1,Nbin], 'Ytick',[1,Nbin]); 
            %xtickangle(30);

            subtightplot(2,3,5, opt{:});  
            imagesc( XCpeak ); axis square;
            colorbar; caxis([0.5,1]);
            set(gca,'Xtick',[1,Nbin], 'Ytick',[1,Nbin]); 
            %xtickangle(30);
            xlabel('Tile'); %ylabel('Tile'); 
            title( sprintf('Peak Cross-Correlation (min = %2.2f)', minXC ));  %title('Peak Cross-Correlation'); 

            subtightplot(2,3,6, opt{:})
            imagesc( XCpeakLag ); hold on;
            %if ~isempty(tempPairInd), plot( tempPair(:,2), tempPair(:,1) ,'k.' ); end
            axis square;
            CB = colorbar; CB.Label.String = 'Delay (s)';
            set(gca,'Xtick',[1,Nbin], 'Ytick',[1,Nbin]); 
            xlabel('Tile'); %ylabel('Tile'); 
            title(sprintf('Cross-Correlation Peak Latency (max = %2.1f)', maxDelay)); %title('Cross-Correlation Peak Latency');
            impixelinfo;
            toc
        end
    end
    clearvars flatVol flatRow flatCol flatInd preFlatVol firstFrame nanScan tempPreFlatScan tempPeriScan periBoutVol rawTrace normTrace filtTrace normDiff
    if ~isempty(saveName)
        fprintf('\nSaving %s', savePath); % , 'normTrace','normDiff','Tperi','onsetWindow','onsetMethod','
        save(savePath, '-v7.3');   
    end
else
    fprintf('\nLoading %s', savePath);
    load(savePath, 'csdWave');  % 'latencyModel',  , 'columnGroup', 'rowGroup', 'planeGroup'
end
toc;

end


%{
lagMat = nan(Nbin,Nbin)
for b = 1:Nbin
    XCpeakLag(1,b)
    lagMat = 
end
%}