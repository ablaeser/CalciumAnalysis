function fluor = GetROIfluor(expt, sbxInfo, ROI, fluor, deform, loco, varargin) % , T
%
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'ROI', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'deform', @isstruct )
addRequired( IP, 'loco', @isstruct )
addOptional( IP, 'axon', [], @isstruct )
%addRequired( IP, 'T', @iscell )
addParameter( IP, 'lp', 0, @isnumeric) % cutoff frequency (Hz) for Gaussian LP filter
addParameter( IP, 'window', 101, @isnumeric) %rolling percentile window
addParameter( IP, 'basePrct', 10, @isnumeric) %rolling percentile value
addParameter( IP, 'deconvolve', false, @islogical) 
addParameter( IP, 'overwrite', false, @islogical)
addParameter( IP, 'reload', true, @islogical)
parse( IP, expt, sbxInfo, ROI, fluor, deform, loco, varargin{:} );  %
%axon = IP.Results.axon;
lpFreq = IP.Results.lp;
windowSize = IP.Results.window;
deconvToggle = IP.Results.deconvolve;
if mod(windowSize,2) ==0,  windowSize = windowSize+1; end
basePrct = IP.Results.basePrct; 
overwrite = IP.Results.overwrite;
reload = IP.Results.reload;
if lpFreq > expt.scanRate, lpFreq = expt.scanRate/(2*pi); end
savePath = strcat( expt.dir, expt.name, '_ROIfluor.mat' );
tic

if ~isfield(fluor(1).F, 'ROI') || overwrite || reload
    if ~exist(savePath, 'file') || overwrite
        Nroi = numel(ROI);
        Nruns = expt.Nruns;
        allROIind = vertcat(ROI.ind); allNPind = vertcat(ROI.neuropil);
        % Get raw signals from ROI and neuropil ROI
        fprintf('\nExtracting traces');
        Fall = nan(expt.totScan, 1); Fback = nan(expt.totScan, 1); Froi = nan(expt.totScan, Nroi); Fnp = nan(expt.totScan, Nroi);
        w = waitbar(0,'Extracting ROI and neuropil traces...'); % w = 
        for s = 1:expt.totScan%-1
            tempVol = double( readSBX(expt.sbx.reg, sbxInfo, s, 1, 1, []) ); 
            %tempVol = double( readSBX(expt.sbx, sbxInfo, expt.Nplane*(s-1)+1, expt.Nplane, 1, []) ); 
            Fall(s) = mean( tempVol(allROIind) ); 
            Fback(s) = mean( tempVol(allNPind) ); 
            for roi = 1:Nroi
                Froi(s,roi) = mean( tempVol(ROI(roi).ind) );
                Fnp(s,roi) = mean( tempVol(ROI(roi).neuropil) );
            end
            waitbar(s/expt.totScan);
        end
        delete(w);
        toc
        clearvars tempVol 
        
        % Low-pass filter the ROI and NP signals
         % {
        if lpFreq > 0
            fprintf('\nLow-pass filtering data (lpFreq = %2.2f)', lpFreq);
            gaussFilt = MakeGaussFilt( 4, 0, 1/(2*pi*lpFreq), expt.scanRate, false );
            Froi_lp = nan(size(Froi)); Fnp_lp = nan(size(Fnp)); Fall_lp = nan(size(Fall)); Fback_lp = nan(size(Fback));  %dFF = nan(size(F)); % Fotot = nan(size(F));  dFFtot = nan(size(F)); filtLength = numel(gaussFilt);
            for run = 1:Nruns
                runScans{run} = expt.scanLims(run)+1:expt.scanLims(run+1);
                for roi = 1:Nroi %find(Nnan>0)
                    Froi_lp(runScans{run},roi) = filtfilt(gaussFilt, 1, Froi(runScans{run},roi));
                    Fnp_lp(runScans{run},roi) = filtfilt(gaussFilt, 1, Fnp(runScans{run},roi));
                end
                Fall_lp(runScans{run}) = filtfilt(gaussFilt, 1, Fall(runScans{run}));
                Fback_lp(runScans{run}) = filtfilt(gaussFilt, 1, Fback(runScans{run}));
            end
            %plot([Fall(:,1), Fall_lp(:,1)])
            %plot([Froi(:,1), Froi_lp(:,1)])
        else
            fprintf('\nLow-pass filtering disabled (lpFreq = %2.2f)', lpFreq);
            for run = 1:Nruns,  runScans{run} = expt.scanLims(run)+1:expt.scanLims(run+1);  end
            Froi_lp = Froi;
            Fnp_lp = Fnp;
            Fall_lp = Fall;
            Fback_lp = Fback;
        end
        
        % Subtract neuropil from ROI 
        F = Froi_lp - Fnp_lp + mean(Fnp_lp,1,'omitnan');  % Froi - Fnp + mean(Fnp,1,'omitnan'); % plot( F(:,1) );
        Ftot = Fall_lp - Fback_lp + mean(Fback_lp,1,'omitnan');
        
        % Calculate baseline signals Fo using rolling percentile
        fprintf('\nCalculating baseline signal');
        Fo = nan(size(F));
        for run = 1:Nruns
            Fo(runScans{run},:) = MovingPercentile(F(runScans{run},:), basePrct, windowSize, 'pre'); % high-pass filter movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
        end
        Fotot = MovingPercentile(Ftot, basePrct, windowSize, 'pre'); 
        
        % Censor scans where each ROI's planes had a bad deformation value
        transAPcat = cat(1, deform.transAP);
        nan_scan = cell(1,Nroi);
        for roi = 1:Nroi
            nan_scan{roi} = find( isnan(sum(transAPcat(:,ROI(roi).zUse),2) ) )';
            F(nan_scan{roi},roi) = NaN;
            fprintf('\nROI %i: censored %i scans with bad deformation values', roi, numel(nan_scan) );
        end
        Nnan = cellfun(@numel, nan_scan);
        
        % Use the censored, neuropil-subtracted, filtered, within-run signal to normalize
        %{
        if lpFreq > 0
            fprintf('\nLow-pass filtering data (lpFreq = %2.2f)', lpFreq);
            gaussFilt = MakeGaussFilt( 4, 0, 1/(2*pi*lpFreq), expt.scanRate, false );
            Ffilt = nan(size(F)); Fo = nan(size(F));  %dFF = nan(size(F)); % Fotot = nan(size(F));  dFFtot = nan(size(F)); filtLength = numel(gaussFilt);
            for run = 1:Nruns
                runScans{run} = expt.scanLims(run)+1:expt.scanLims(run+1);
                for r = 1:Nroi %find(Nnan>0)
                    %nanCC = bwconncomp( isnan( F(runScans{run},r) ) );
                    tempCC = bwconncomp( ~isnan( F(runScans{run},r)) );
                    %CClength = cellfun(@numel, tempCC.PixelIdxList);
                    for c = 1:tempCC.NumObjects
                        %plot(F(tempConnComp.PixelIdxList{c},r)); hold on; plot(filtfilt(gaussFilt, 1, F(tempConnComp.PixelIdxList{c},r)));
                        try
                            Ffilt(expt.scanLims(run)+tempCC.PixelIdxList{c},r) = filtfilt(gaussFilt, 1, F(expt.scanLims(run)+tempCC.PixelIdxList{c},r)); % filtfilt can't handle NaN
                        catch
                            fprintf('\n[run, roi, stretch] = [%i, %i, %i]: gaussian filtering failed!', run, r, c)
                        end
                    end
                end
                Fo(runScans{run},:) = MovingPercentile(Ffilt(runScans{run},:), basePrct, windowSize, 'pre'); % high-pass filter movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
            end
            %Fo = MovingPercentile(F, basePrct, windowSize, 'pre'); % high-pass filter movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
        else
            fprintf('\nLow-pass filtering disabled (lpFreq = %2.2f)', lpFreq);
            Ffilt = F;
            Fo = nan(size(F));  %dFF = nan(size(F)); % Fotot = nan(size(F));  dFFtot = nan(size(F)); filtLength = numel(gaussFilt);
            for run = 1:Nruns
                runScans{run} = expt.scanLims(run)+1:expt.scanLims(run+1);
                Fo(runScans{run},:) = MovingPercentile(Ffilt(runScans{run},:), basePrct, windowSize, 'pre'); % high-pass filter movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
            end
        end
        %}
        dFF = (F-Fo)./Fo; % (Ffilt-Fo)./Fo; %dFF(tempCC.PixelIdxList{c},:) = (F-Fo)./Fo;
        dFFtot = (Ftot-Fotot)./Fotot;

        % Deconvolve dF/F
        if deconvToggle % && expt.Nplane == 1
            fprintf('\nDeconvolving dF/F')
            tic
            runAct = cell(1,Nruns); %runScans = cell(1,Nruns);
            parfor run = 1:Nruns % 
                %fprintf('run = %i / %i', run, Nruns) 
                %runScans{run} = expt.scanLims(run)+1:expt.scanLims(run+1);
                %w = waitbar(0, 'Deconvolving censored dF/F'); %fprintf('\nDeconvolving censored dF/F'); sprintf('Deconvolving censored dF/F (run %i of %i)', run, Nruns)
                runAct{run} = nan(numel(runScans{run}), Nroi);
                for roi = 1:Nroi %for r = % flip(1:Nroi) % 
                    try
                        goodRunFrames = find(~isnan(dFF(runScans{run},roi)))'; %intersect(6, );
                        dFFrunCens = dFF(runScans{run}(goodRunFrames),roi);
                        [~, actRunCens, ~] = deconvolveCa( dFFrunCens, 'exp2', 'foopsi', 'optimize_pars', 'optimize_smin', 'maxIter',20, 'window',300 );  %ar1  deconvOpts(r)
                        runAct{run}(goodRunFrames,roi) = actRunCens; % 
                        %activity(goodFrames,r) = actRunCens;
                        %figure; yyaxis left; plot( dFFrunCens ); yyaxis right; plot( actRunCens ); pause;
                    catch
                        fprintf('\nROI %i:  Deconvolution failed', roi);
                        %rFailed = [rFailed, r];
                    end
                    %waitbar(r/Nroi, w);
                    %delete(w);
                end
                %toc
            end
            activity = vertcat(runAct{:});
            activity(activity<0) = 0;
            %{
            w = waitbar(0, 'Deconvolving censored dF/F');
            %fprintf('\nDeconvolving censored dF/F');
            for r = flip(1:Nroi)
                try
                    goodFrames = find(~isnan(dFF(:,r)))';
                    dFFcens = dFF(goodFrames,r);
                    [~, actCens, ~] = deconvolveCa( dFFcens, 'exp2', 'foopsi', 'optimize_pars', 'optimize_smin', 'maxIter',20, 'window',300 );  %ar1  deconvOpts(r)
                    activity(goodFrames,r) = actCens;
                    %figure; yyaxis left; plot( dFF(:,r) ); yyaxis right; plot( activity(:,r) ); 
                catch
                    fprintf('\nROI %i:  Deconvolution failed', r);
                    rFailed = [rFailed, r];
                end
                waitbar(r/Nroi, w);
            end
            delete(w);
            %}
            toc
        else
            fprintf('\nDeconvolution disabled')
            activity = nan(size(dFF,1),size(dFF,2));
        end
        %figure; yyaxis left; plot( dFF(:,1) ); yyaxis right; plot( activity(:,1) ); pause;
        %}

        %{
        % Find ROIs' brightest scans and use them to make a projection image
        figure('WindowState','maximized');
        Npeaks = 30;
        for r = 1:Nroi
            [~,brightScans] = sort( dFF(:,r), 'descend', 'MissingPlacement','last'); brightScans = brightScans(1:Npeaks);
            for p = flip(1:Npeaks)
                peakVol(:,:,:,p) = pipe.io.sbxRead(expt.sbx, expt.Nplane*(brightScans(p)-1)+1, expt.Nplane, 1, []);
            end
            meanPeakVol = mean(peakVol, 4);
            ROI(r).peakProj = max( meanPeakVol(:,:,[ceil(ROI(r).box_z(1)), floor( ROI(r).box_z(2) )]), [], 3); 
            %cropProj = max(meanPeakVol(ROI(r).box_y(1):ROI(r).box_y(2), ROI(r).box_x(1):ROI(r).box_x(2), [ceil(ROI(r).sub(1).box_z(1)), floor( ROI(r).sub(1).box_z(2) )]), [], 3); 
            %{
            cropPrctile = prctile( cropProj(:), [10, 99.5] );
            subplot(2,1,1); cla;
            plot( dFF(:,r) ); hold on;
            plot( brightScans, dFF(brightScans,r), 'ko' ); hold off;
            axis tight;
            subplot(2,1,2); cla;
            imshow( cropProj, cropPrctile ); hold on; % 
            for s = 1:ROI(r).Nsub
                plot( ROI(r).edge(:,2)-ROI(r).box_x(1)+1, ROI(r).edge(:,1)-ROI(r).box_y(1)+1, '.', 'MarkerSize',2 ); %, 'Color',tempColors(s,:)
            end
            impixelinfo;
            pause;
            %}
        end
        %}
        fprintf('\nSaving %s', savePath);
        save(savePath, 'Froi','Fnp','Froi_lp','Fnp_lp','F','Fo','dFF','activity','Fall','Fback','Fall_lp','Fback_lp','Ftot','Fotot','dFFtot','nan_scan','Nroi','Nruns','lpFreq', '-v7.3'); % ,'runScans' ,'rFailed'
    else
        fprintf('\nLoading %s', savePath);
        load(savePath, 'Froi','Fnp','Froi_lp','Fnp_lp','F','Fo','dFF','activity','Fall','Fback','Fall_lp','Fback_lp','Ftot','Fotot','dFFtot','nan_scan','Nroi','Nruns','lpFreq' ); % ,'nan_scan'
        % 'Fall','Fback','Froi','Fnp','Ftot','F','Fo','Fotot','dFF','dFFtot','activity','nan_scan','Nroi','Nruns','lpFreq'
    end

    % Break fluor signals up by run and normalize using periods of stillness
    [fluor, stillScans] = StillNormalizedFluor(expt, fluor, loco, Fall, Fback, Ftot, Fotot, dFFtot, Froi, Fnp, F, Fo, dFF, activity);
    toc
else
    fprintf('\nROI fluor already found!');
end

%{
% Break fluor signals up by run
    fprintf('\nBreaking fluor signals up by run...  ')
    scanLims = expt.scanLims;
    for runs = 1:expt.Nruns
        % Pooled ROI
        fluor(runs).Froi.all = Fall(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).Fnp.all = Fback(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).F.all = Ftot(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Fo.all = Fotot(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).dFF.all = dFFtot(scanLims(runs)+1:scanLims(runs+1),:);      
        fluor(runs).z.all = (fluor(runs).dFF.all - median(fluor(runs).dFF.all(loco(runs).stateDown == 1,:), 1, 'omitnan'))./std(fluor(runs).dFF.all(loco(runs).stateDown == 1,:), 1, 'omitnan');% normalize(fluor(runs).dFF.all, 1);
        % Individual ROI
        fluor(runs).Froi.ROI = Froi(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Froi.med = median(fluor(runs).Froi.ROI, 1, 'omitnan'); %(loco(runs).stateDown == 1,:)
        fluor(runs).Fnp.ROI = Fnp(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Fnp.med = median(fluor(runs).Fnp.ROI, 1, 'omitnan'); %(loco(runs).stateDown == 1,:)
        fluor(runs).Froi.SNR = (fluor(runs).Froi.med - fluor(runs).Fnp.med)./fluor(runs).Fnp.med;
        fluor(runs).F.ROI = F(scanLims(runs)+1:scanLims(runs+1),:);        
        fluor(runs).Fo.ROI = Fo(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).dFF.ROI = dFF(scanLims(runs)+1:scanLims(runs+1),:);
        %{
        testSignal = fluor(runs).dFF.ROI(:,1);
        std( testSignal(), 'omitnan' )
        
        sp(1) = subplot(3,1,1); plot(testSignal);
        sp(2) = subplot(3,1,2); plot(normalize(testSignal,1)); hold on; %plot(testSignal - mean(testSignal,'omitnan')/std(testSignal,'omitnan'));
        line([0,size(testSignal,1)], [0,0],'color','k');
        
        stillMean = median(testSignal(loco(runs).stateDown == 1), 'omitnan'); %mean(testSignal(loco(runs).stateDown == 1), 'omitnan');
        stillStd = std(testSignal(loco(runs).stateDown == 1), 'omitnan');
        sp(3) = subplot(3,1,3); plot((testSignal-stillMean)/stillStd); hold on; 
      
        line([0,size(testSignal,1)], [0,0],'color','k');
        linkaxes(sp,'xy');
        axis tight;
        %}
        fluor(runs).z.ROI = (fluor(runs).dFF.ROI - median(fluor(runs).dFF.ROI(loco(runs).stateDown == 1,:), 1, 'omitnan'))./std(fluor(runs).dFF.ROI(loco(runs).stateDown == 1,:), 1, 'omitnan'); %normalize(fluor(runs).dFF.ROI, 1); % 
        fluor(runs).act.ROI = activity(scanLims(runs)+1:scanLims(runs+1),:);
        % Correlation between ROIs for this run
        fluor(runs).Froi.corr = corr( fluor(runs).Froi.ROI, 'Rows','complete' );
        fluor(runs).Fnp.corr = corr( fluor(runs).Fnp.ROI, 'Rows','complete' );
        fluor(runs).F.corr = corr( fluor(runs).F.ROI, 'Rows','complete' );
        fluor(runs).Fo.corr = corr( fluor(runs).Fo.ROI, 'Rows','complete' );
        fluor(runs).dFF.corr = corr( fluor(runs).dFF.ROI, 'Rows','complete' );
        fluor(runs).act.corr = corr( fluor(runs).act.ROI, 'Rows','complete' );
        %{
        for r = 1:Nroi
            subplot(3,2,[1,3,5]);
            imshow( ROI(r).mask_2d, [] ); hold on;
            title( sprintf('ROI %i',r) );
            impixelinfo
            sp(1) = subplot(3,2,2); cla;
            plot( fluor(runs).F.ROI(:,r) ); hold on;
            plot( fluor(runs).Fo.ROI(:,r) ); ylabel('F');
            title(sprintf('[run, ROI] = [%i, %i]', runs, r));
            sp(2) = subplot(3,2,4); cla;
            plot( fluor(runs).dFF.ROI(:,r) ); ylabel('dF/F');
            sp(3) = subplot(3,2,6); cla;
            plot( fluor(runs).act.ROI(:,r) ); ylabel('Activity');
            linkaxes(sp,'x'); xlim([-Inf,Inf]);
            pause;
        end
        %}
    end
%}