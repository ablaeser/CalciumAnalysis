function csdWave = GetNeuropilWaveSpeed(expt, catInfo, periBout, voxInd, varargin)
% Track how fast a CSD wave washes across a set of specific volumes (ROIs, axons or neuropil)
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
addRequired( IP, 'periBout', @isstruct )
addRequired( IP, 'voxInd', @iscell )
addParameter( IP, 'Tshift', 0, @isnumeric )
addParameter( IP, 'window', [-6,14], @isnumeric) % time window (seconds) in which to detect onset of wave
addParameter( IP, 'bin', 40, @isnumeric) % width of binning (pixels)
%addParameter( IP, 'minSep', 200/expt.scanRate, @isnumeric ) % 150
%addParameter( IP, 'minDFF', 0.1, @isnumeric )
%addParameter( IP, 'minXC', 0.5, @isnumeric )
%addParameter( IP, 'maxDelay', 15, @isnumeric )
%addParameter( IP, 'method', 'peak', @isnumeric )
addParameter( IP, 'lp', 2, @isnumeric) % cutoff frequency (Hz) for Gaussian LP filter
addParameter( IP, 'minR2', 0.5, @isnumeric )
addParameter( IP, 'maxTau', 2, @isnumeric )
addParameter( IP, 'show', false, @islogical)
addParameter( IP, 'save', '', @ischar)
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, catInfo, periBout, voxInd, varargin{:} );  %, fluor
binSize = IP.Results.bin;
csdWave.window = IP.Results.window;
%minDFF = IP.Results.minDFF;
%minXC = IP.Results.minXC;
%minSep = IP.Results.minSep;
%maxDelay = IP.Results.maxDelay;
%delayMethod = IP.Results.method;
lpFreq = IP.Results.lp;
minR2 = IP.Results.minR2;
maxTau = IP.Results.maxTau;
csdWave.Tshift = IP.Results.Tshift;
show = IP.Results.show;
saveName = IP.Results.save;
overwrite = IP.Results.overwrite;

if ~isempty(saveName),  savePath = sprintf('%s%s_wave_%s.mat', expt.dir, expt.name, saveName);  end
if isempty(saveName) || ~exist(savePath,'file') || overwrite
    if lpFreq > expt.scanRate, lpFreq = expt.scanRate/(2*pi); end
    gaussSigma = 1/(2*pi*lpFreq);
    gaussFilt = MakeGaussFilt( 4, 0, gaussSigma, expt.scanRate, false ); %  gaussFilt = MakeGaussFilt( 4, 0, 1/expt.scanRate, expt.scanRate, true );
    modelSigmoid = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4) ); % beta = [factor, growth rate, threshold, constant response offset]
    % Determine initial parameters for sigmoidal fit
    betaInit(1) = 1; % Saturation response
    betaInit(2) = 1; % time constant
    betaInit(3) = 0; % Onset time
    betaInit(4) = 0; % Constant offset response
    for b = 1:periBout.Nbout
        % Get peri-bout volume scans
        periBoutScans = sum(expt.Nscan(1:expt.csd-1)) + periBout.scan{b};
        periBoutVol = double(readSBX(expt.sbx, catInfo, periBoutScans(1), periBout.Nscan(b), 1, []));
        %periBoutVol = double(readSBX(expt.sbx, catInfo, expt.Nplane*(periBoutScans(1)-1)+1, numel(periBoutScans)*expt.Nplane, 1, [])); % double(pipe.io.sbxRead(expt.sbx, expt.Nplane*(periCSDscans(1)-1)+1, numel(periCSDscans)*expt.Nplane, 1, [])); %zeros( )
        %periBoutVol = reshape( periBoutVol, expt.Nrow, expt.Ncol, expt.Nplane, periBout.Nscan(b) );

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
        %saveastiff( uint16(flatVol), sprintf('%s%s_flatVol.tif', expt.dir, expt.name ) )
        %flatVol(isnan(flatVol)) = 0;
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
        
        % Find bad scans
        csdWave.badScan = find(isnan(sum(periBout.fluor{1},2)))';
        csdWave.Nbad = numel(csdWave.badScan);
        if csdWave.Nbad > 0, warning('Found %i bady registered scans - excluded from analysis', csdWave.Nbad); end

        % Get periBout timing and define window (if not already defined)
        Tperi = periBout.T{b} - periBout.Tstart(b) - csdWave.Tshift; 
        csdWave.scans = setdiff( find(Tperi >= csdWave.window(1) & Tperi <= csdWave.window(2)), csdWave.badScan );
        csdWave.T = Tperi(csdWave.scans);
        csdWave.preScan = find(csdWave.T < 0)'; % 1:floor(expt.scanRate*preTime); %
        csdWave.postScan = find(csdWave.T >= 0)'; %ceil(expt.scanRate*preTime):ceil(expt.scanRate*preTime)+floor(expt.scanRate*postTime); % 
        
        % Bin traces into x-y boxes
        tic
        fprintf('\nBinning pixels (binSize = %i)', binSize)
        minX = min(traceXYZ{b}(:,1));  maxX = max(traceXYZ{b}(:,1));  minY = min(traceXYZ{b}(:,2));  maxY = max(traceXYZ{b}(:,2));
        csdWave.binSize = binSize;
        csdWave.xBinLim = minX:binSize:maxX;
        csdWave.yBinLim = minY:binSize:maxY;
        csdWave.NbinX = numel(csdWave.xBinLim)-1; csdWave.NbinY = numel(csdWave.yBinLim)-1;
        csdWave.Nbin = csdWave.NbinX*csdWave.NbinY; %size(binTrace,2);
        csdWave.binTrace = nan(csdWave.NbinX, csdWave.NbinY, periBout.Nscan(b)); csdWave.binCent = nan(csdWave.NbinX, csdWave.NbinY, 3);
        for x = 1:csdWave.NbinX
            for y = 1:csdWave.NbinY
                boxTraceInd = find( traceXYZ{b}(:,1) >= csdWave.xBinLim(x) & traceXYZ{b}(:,1) < csdWave.xBinLim(x+1) & traceXYZ{b}(:,2) >= csdWave.yBinLim(y) & traceXYZ{b}(:,2) < csdWave.yBinLim(y+1) );
                csdWave.binTrace(x,y,:) = filtfilt(gaussFilt, 1, mean(rawTrace{b}(:,boxTraceInd), 2, 'omitnan') ); % filtTrace
                csdWave.binCent(x,y,:) = expt.umPerPixel*mean(traceXYZ{b}(boxTraceInd,:));
            end
        end       
        csdWave.binTrace = csdWave.binTrace(:,:,csdWave.scans);

        % figure; subplot(2,1,1); imagesc( csdWave.binTraceDFF' ); colorbar; impixelinfo; axis image; subplot(2,1,2); ecdf(csdWave.binTraceDFF(:)); hold on; line(minDFF*[1,1], [0,1], 'color','k'); xlabel('dF/Fo'); 
        csdWave.normTrace = (csdWave.binTrace - median( csdWave.binTrace(:,:,csdWave.preScan), 3 )); % min(csdWave.binTrace(:,:,csdWave.preScan),[],3))
        csdWave.normTrace = csdWave.normTrace./max(csdWave.normTrace(:,:,csdWave.postScan),[],3);
        %{
        %[xBad, yBad] = ind2sub( size(csdWave.normTrace), find( max( csdWave.normTrace, [], 3 ) > 1 | min( csdWave.normTrace, [], 3 ) < 0 | csdWave.binTraceDFF < minDFF ) ); %[xBad, yBad]
        %Nbad = numel(xBad);
        %fprintf('\n%i of %i  bins excluded from analysis', Nbad, csdWave.Nbin);
        %min( csdWave.normTrace, [], 3 )
        
        %figure;
        for p = 1:numel(xBad)
            subplot(1,2,1); plot( squeeze(csdWave.binTrace(xBad(p),yBad(p),:))); title( sprintf('dF/F = %2.2f', csdWave.binTraceDFF(xBad(p),yBad(p))) );
            subplot(1,2,2); plot(squeeze(csdWave.normTrace(xBad(p),yBad(p),:))); pause;
            csdWave.normTrace(xBad(p),yBad(p),:) = 0;
        end
        %}
        %saveastiff( uint16(permute(100*csdWave.normTrace, [2,1,3])), sprintf('%s%s_csdWave.normTrace.tif', expt.dir, expt.name ) )

        % Estimate wave onset by fitting to sigmoid            
        %close all; figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
        csdWave.binR2 = nan(csdWave.NbinX, csdWave.NbinY); csdWave.binOnsetTime = nan(csdWave.NbinX, csdWave.NbinY); csdWave.binTimeConst = nan(csdWave.NbinX, csdWave.NbinY); 
        csdWave.sigmoidData = nan(csdWave.Nbin, 5); k =0;
        for x = 1:csdWave.NbinX
            for y = 1:csdWave.NbinY
                sigmFit = fitnlm( csdWave.T, squeeze(csdWave.normTrace(x,y,:)), modelSigmoid, betaInit );
                fitResult.R2 =  sigmFit.Rsquared.Adjusted; % goodness of fit
                fitResult.Tonset = sigmFit.Coefficients.Estimate(3);
                fitResult.pTonset = sigmFit.Coefficients.pValue(3);
                fitResult.tau = sigmFit.Coefficients.Estimate(2);
                fitResult.pTau = sigmFit.Coefficients.pValue(2);
                fitResult.factor = sigmFit.Coefficients.Estimate(1);
                fitResult.pFactor = sigmFit.Coefficients.pValue(1);
                fitResult.offset = sigmFit.Coefficients.Estimate(4);
                fitResult.pOffset = sigmFit.Coefficients.pValue(4);
                csdWave.binR2(x,y) = fitResult.R2;
                if fitResult.R2 > minR2 && fitResult.tau < maxTau && fitResult.tau >= 0
                    csdWave.binTimeConst(x,y) = fitResult.tau;
                    csdWave.binOnsetTime(x,y) = fitResult.Tonset; 
                    k = k+1;
                    csdWave.sigmoidData(k,:) = [csdWave.binCent(x,y,1), csdWave.binCent(x,y,2), fitResult.Tonset, fitResult.tau, fitResult.R2];
                end
                %{
                plot( csdWave.T, squeeze(csdWave.normTrace(x,y,:)) ); hold on;
                plot( csdWave.T, predict(sigmFit, csdWave.T), 'r' );
                title(sprintf('[x,y] = [%i, %i] R^2 = %2.2f, onset = %2.1f  (p=%2.4f),  time constant = %2.3f (p=%2.4f)', ...
                    x,y, fitResult.R2, fitResult.Tonset, fitResult.pTonset, fitResult.tau, fitResult.pTau) );
                pause; cla;
                %}
            end
        end
        csdWave.sigmoidData(k+1:end,:) = [];
        
        % Calculate speed
        csdWave.speedMdl = fitlm( sqrt( sum(csdWave.sigmoidData(:,1:2).^2, 2) ), csdWave.sigmoidData(:,3), 'RobustOpts',true);
        csdWave.speed = abs(1/csdWave.speedMdl.Coefficients.Estimate(2));
        % CALCULATING VELOCITY THIS WAY IS SUSPECT, ESPECIALLY THE Y COMPONENT
        csdWave.velocityMdl = fitlm( csdWave.sigmoidData(:,1:2), csdWave.sigmoidData(:,3), 'RobustOpts',true);
        csdWave.velocity = [0,0]; %
        csdWave.pVelocity(1) = csdWave.velocityMdl.Coefficients.pValue(2);  
        csdWave.pVelocity(2) = csdWave.velocityMdl.Coefficients.pValue(3); 
        csdWave.velocity(1) = 1./csdWave.velocityMdl.Coefficients.Estimate(2); %if csdWave.pVelocity(1) < 0.01,  end
        %if csdWave.pVelocity(2) < 0.01 %&& abs(1./csdWave.velocityMdl.Coefficients.Estimate(2)) <= abs(csdWave.speed) %end  
        csdWave.velocity(2) = -1./csdWave.velocityMdl.Coefficients.Estimate(3);  % y axis is flipped

        % Use contour to get angle of wavefront at median timepoint
        csdWave.contour = contourc( csdWave.binOnsetTime', median(csdWave.binOnsetTime(:), 'omitnan')*[1,1] ); % flip( ,1)
        frontNodes = find( csdWave.contour(1,:) == csdWave.contour(1,1) ); % prevent gaps in the front 
        if numel(frontNodes) > 1
            [~,maxNodeInd] = max(csdWave.contour(2,frontNodes));
            try
                csdWave.contour = csdWave.contour(:,frontNodes(maxNodeInd):frontNodes(maxNodeInd+1)-1);
            catch
                csdWave.contour = csdWave.contour(:,frontNodes(maxNodeInd):end);
            end
        end
        if range(csdWave.contour(1,2:end)) >= range(csdWave.contour(2,2:end))
            csdWave.frontFit = fitlm( csdWave.contour(1,2:end), csdWave.contour(2,2:end) ); % , 'RobustOpts',true
            frontEnds = [csdWave.contour(1,[2,end])', predict(csdWave.frontFit, csdWave.contour(1,[2,end])')]; % %[x1, y1; x2, y2]
        else
            csdWave.frontFit = fitlm( csdWave.contour(2,2:end), csdWave.contour(1,2:end), 'RobustOpts',true );
            frontEnds = [predict(csdWave.frontFit, csdWave.contour(2,[2,end])'), csdWave.contour(2,[2,end])']; %[x1, y1; x2, y2]
        end
        csdWave.frontVec = diff(frontEnds, 1, 1);
        csdWave.propVec = [sign(csdWave.velocity(1))*abs(csdWave.frontVec(2)), sign(csdWave.velocity(2))*abs(csdWave.frontVec(1))];  % direction of propagation is perpendicular to wavefront AND Y axis is reversed
        csdWave.angle = cart2pol( csdWave.propVec(1), csdWave.propVec(2) ); 

        %{
        figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
        subplot(2,3,1); imagesc(csdWave.binR2' ); title('R^2') ; colorbar; axis image;
        subplot(2,3,4); ecdf(csdWave.binR2(:)); xlabel('R^2');
        subplot(2,3,2); imagesc(csdWave.binOnsetTime'); title('Onset Time'); colorbar; axis image; hold on; impixelinfo;
        contour(csdWave.binOnsetTime')
        %plot( csdWave.contour(1,2:end), csdWave.contour(2,2:end), 'k' );
        %line( csdWave.contour(1,[2,end]), predict(csdWave.frontFit, csdWave.contour(1,[2,end])'), 'linestyle','--' )
        subplot(2,3,5); ecdf(csdWave.binOnsetTime(:)); xlabel('Onset Time (s)');
        subplot(2,3,3); imagesc(csdWave.binTimeConst' ); title('Time Constant') ; colorbar; axis image; impixelinfo
        subplot(2,3,6); ecdf(csdWave.binTimeConst(:)); xlabel('Time Constant (s)');       
        %}
    end
    close all;
    clearvars flatVol flatRow flatCol flatInd preFlatVol firstFrame nanScan tempPreFlatScan tempPeriScan periBoutVol rawTrace normTrace filtTrace normDiff
    if ~isempty(saveName)
        fprintf('\nSaving %s\n', savePath); 
        save(savePath, '-v7.3');   
    end
else
    fprintf('\nLoading %s\n', savePath);
    load(savePath, 'csdWave');  
end

if show
    binTraceXsort = nan( size(csdWave.normTrace,3), csdWave.Nbin);  binTraceYsort = nan( size(csdWave.normTrace,3), csdWave.Nbin);
    k = 0;
    for x = 1:csdWave.NbinX
        for y = 1:csdWave.NbinY
            if sum(csdWave.normTrace(x,y,:),3) ~= 0 
                k = k+1;
                binTraceXsort(:,k) = squeeze(csdWave.normTrace(x,y,:)); 
            end
        end
    end
    binTraceXsort(:,isnan(sum(binTraceXsort,1))) = [];
    k = 0;
    for y = 1:csdWave.NbinY
        for x = 1:csdWave.NbinX
            if sum(csdWave.normTrace(x,y,:),3) ~= 0
                k = k+1;
                binTraceYsort(:,k) = squeeze(csdWave.normTrace(x,y,:)); 
            end
        end
    end
    binTraceYsort(:,isnan(sum(binTraceYsort,1))) = [];

    opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
    FS = 10;
    close all; clearvars sp;  
    XtickZero = find(csdWave.T == 0);
    XtickWindow = 1:XtickZero-1:numel(csdWave.T); 

    figure('Units','normalized', 'OuterPosition',[0,0,1,1]);
    % Sorted tile data heatmaps
    subtightplot(2,4,1,opt{:}); 
    imagesc( binTraceXsort' ); %axis image;
    cb = colorbar; cb.Label.String = 'Normalized'; caxis([0,1]);
    ylabel('Tile'); title(sprintf('Binned data (width = %i pixels)', binSize));  %xlabel('Peri-CSD Time (s)'); 
    set(gca, 'Ytick',[1,size(binTraceXsort,2)], 'YtickLabel',{'Leftmost','Rightmost'}, 'Xtick',XtickWindow, 'XtickLabel',sprintfc('%2.1f', csdWave.T(XtickWindow)), 'TickDir','out' ); % 

    subtightplot(2,4,5,opt{:}); 
    imagesc( binTraceYsort' ); %axis square;%axis image;
    cb = colorbar; cb.Label.String = 'Normalized'; caxis([0,1]);
    xlabel('Peri-CSD Time (s)');  ylabel('Tile');
    set(gca, 'Ytick',[1,size(binTraceYsort,2)], 'YtickLabel',{'Topmost','Bottommost'}, 'Xtick',XtickWindow, 'XtickLabel',sprintfc('%2.1f', csdWave.T(XtickWindow)), 'TickDir','out' ); 
    xlim([-Inf,Inf]);

    % Sigmoidal fit results
    subtightplot(3,4,2,opt{:}); 
    imagesc(csdWave.binOnsetTime'); title('Binned Map: Onset Time (sigmoidal fit)') ; colorbar; axis image; hold on;     
    plot( csdWave.contour(1,2:end), csdWave.contour(2,2:end), 'k' );
    line( [csdWave.contour(1,2), csdWave.contour(1,2)+csdWave.frontVec(1)], [csdWave.contour(2,2), csdWave.contour(2,2)+csdWave.frontVec(2)], 'linestyle','--')
    set(gca,'Xtick',[], 'Ytick',[]);

    subtightplot(3,4,6,opt{:}); 
    imagesc(csdWave.binTimeConst'); title('Time Constant') ; colorbar; axis image;
    set(gca,'Xtick',[], 'Ytick',[]);

    subtightplot(3,4,10,opt{:}); 
    imagesc(csdWave.binR2'); title('R^2') ; colorbar; axis image;
    set(gca,'Xtick',[], 'Ytick',[]);

    % Linear fits 
    subtightplot(2,4,4,opt{:});
    plot( sqrt( sum(csdWave.sigmoidData(:,1:2).^2, 2) ), csdWave.sigmoidData(:,3), '.' ); hold on;
    plot( sqrt( sum(csdWave.sigmoidData(:,1:2).^2, 2) ), predict(csdWave.speedMdl, sqrt( sum(csdWave.sigmoidData(:,1:2).^2, 2) )), 'r' )
    xlabel('Radial Distance (um)'); ylabel('Onset time (s)');
    axis square;

    subtightplot(2,4,3,opt{:}); 
    plot( csdWave.sigmoidData(:,1), csdWave.sigmoidData(:,3), '.' ); hold on;
    plot( csdWave.sigmoidData(:,1), predict(csdWave.velocityMdl, [csdWave.sigmoidData(:,1),zeros(size(csdWave.sigmoidData,1),1)]), 'r' )
    xlabel('X Position (um)'); ylabel('Onset time (s)');
    axis square;

    subtightplot(2,4,7,opt{:}); 
    plot( csdWave.sigmoidData(:,2), csdWave.sigmoidData(:,3), '.' ); hold on;
    plot( csdWave.sigmoidData(:,2), predict(csdWave.velocityMdl, [zeros(size(csdWave.sigmoidData,1),1), csdWave.sigmoidData(:,2)]), 'r' )           
    xlabel('Y Position (um)'); ylabel('Onset time (s)');
    axis square;

    % Polar plot of wave velocity
    subtightplot(2,4,8,opt{:});
    polarplot( csdWave.angle*[1,1], csdWave.speed*[0,1], 'LineWidth',3, 'Color','r' ); hold off;
    %axis square;
    title('Estimated CSD Wave Velocity (um/s)', 'FontSize',FS);          
    impixelinfo;
end

toc;
end