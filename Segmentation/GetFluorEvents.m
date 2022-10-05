function [fluorEvent, Nevent, actThresh] = GetFluorEvents(expt, T, fluor, varargin)
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @iscell )
addRequired( IP, 'fluor', @isstruct ) 
addParameter( IP, 'peri', 60, @isnumeric )
addParameter( IP, 'zThresh', 2, @isnumeric )
addParameter( IP, 'minSep', 6, @isnumeric ) % seconds
addParameter( IP, 'show', false, @islogical )
parse( IP, expt, T, fluor, varargin{:} );  
%periTime = IP.Results.peri;
%shortTime = 2; %longTime = 5;
zThresh = IP.Results.zThresh;
minSep = IP.Results.minSep;
show = IP.Results.show;

% Detect high sigma events for each ROI
Nrun = numel(fluor);
%periScan = round(periTime*expt.scanRate);
minSepScan = round(minSep*expt.scanRate);
%shortScan = round(shortTime*expt.scanRate);
%longScan = round(longTime*expt.scanRate);
Nroi = expt.Nroi;
Nevent = nan(Nrun, Nroi);
NscanRun = expt.Nscan;
fluorEvent = cell(Nrun,Nroi); actThresh = cell(Nrun,1);
tic
if ~all(isnan(sum(fluor(1).act.ROI, 1, 'omitnan'))) % use activity if available
    threshPrctile = normcdf(zThresh);
    parfor run = 1:Nrun
        % Determine activity thresholds
        [actLevel,cumFreq] = ecdfCol(fluor(run).act.ROI); % , 'Activity'
        for roi = 1:Nroi
            [~,minInd] = min( abs(cumFreq(:,roi) - threshPrctile), [], 1 );
            actThresh{run}(roi) = actLevel(minInd,roi); %cumFreq(minInd)
        end
        %actThresh{run} = prctile(fluor(run).act.ROI, threshPrctile, 1);
        
        for roi = 1:Nroi
            [evtActPeaks, evtPeakScan] = findpeaks(fluor(run).act.ROI(:,roi), 'MinPeakHeight', actThresh{run}(roi), 'MinPeakDistance',minSepScan, 'SortStr','descend');
            Nevent(run,roi) = numel(evtActPeaks);
            for evt = 1:Nevent(run,roi)
                tempScan = evtPeakScan(evt)-minSepScan:evtPeakScan(evt)+minSepScan;
                tempScan(tempScan > NscanRun(run) | tempScan < 1) = [];
                fluorEvent{run,roi}(evt).periScan = tempScan;
                fluorEvent{run,roi}(evt).Tpeak = T{run}(fluorEvent{run,roi}(evt).periScan);
                fluorEvent{run,roi}(evt).Tperi = T{run}(fluorEvent{run,roi}(evt).periScan) - T{run}(evtPeakScan(evt));
                fluorEvent{run,roi}(evt).z = fluor(run).z.ROI(fluorEvent{run,roi}(evt).periScan, roi);
                fluorEvent{run,roi}(evt).zMax = max( fluorEvent{run,roi}(evt).z );
                fluorEvent{run,roi}(evt).dFF = fluor(run).dFF.ROI(fluorEvent{run,roi}(evt).periScan, roi);
                fluorEvent{run,roi}(evt).dFFmax = max( fluorEvent{run,roi}(evt).dFF );
                fluorEvent{run,roi}(evt).dFFbase = prctile( fluorEvent{run,roi}(evt).dFF, 10); %max( fluor(r).dFF.ROI(fluorEvent{r,roi}(evt).periScan, roi) );
                fluorEvent{run,roi}(evt).act = fluor(run).act.ROI(fluorEvent{run,roi}(evt).periScan, roi);
                fluorEvent{run,roi}(evt).actMax = evtActPeaks(evt); %max( fluorEvent{run,roi}(evt).dFF );
                fluorEvent{run,roi}(evt).actBase = prctile( fluorEvent{run,roi}(evt).act, 10); %max( fluor(r).dFF.ROI(fluorEvent{r,roi}(evt).periScan, roi) );
            end
        end
    end
    toc
end
actThresh = vertcat(actThresh{:});

if show
    close all;
    opt = {[0.08,0.08], [0.08,0.06], [0.08,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
    figure('Units','normalized', 'OuterPosition',[0,0,1,1]); 
    sp(2) = subtightplot(2,1,2,opt{:}); sp(1) = subtightplot(2,1,1,opt{:}); 
    linkaxes(sp,'x');
    for run = 1:Nrun
        for roi = 1:Nroi
            for evt = 1:Nevent(run,roi)
                subplot(sp(1)); cla;
                plot(fluorEvent{run,roi}(evt).Tperi, fluor(run).dFF.ROI(fluorEvent{run,roi}(evt).periScan,roi) ); ylim([-0.2,1]); 
                ylabel('dF/Fo'); xlabel('Absolute Time (s)');
                subplot(sp(2)); cla;
                plot(fluorEvent{run,roi}(evt).Tperi, fluor(run).act.ROI(fluorEvent{run,roi}(evt).periScan,roi) ); hold on;
                line(fluorEvent{run,roi}(evt).Tperi([1,end]), actThresh(run,roi)*[1,1], 'lineStyle','--', 'color','r');
                ylabel('Deconvolved Activity'); xlabel('Peri-Peak Time (s)');
                title( sprintf('[run, roi, event] = [%i, %i, %i]', run, roi, evt) );
                xlim([-Inf,Inf]);
                pause;
            end
        end
    end
end
end

%{
    tempCrossings = fluor(run).z.ROI > zThresh;
    tempEvents = imdilate(tempCrossings, strel('rectangle',[periScan,1]) );

    if show
        clf;
        sp(1) = subplot(4,1,1); imagesc( fluor(run).z.ROI' ); 
        tempPos = get(gca,'Position');
        caxis([-1,3]); 
        colorbar; 
        set(gca,'Position',tempPos);
        sp(2) = subplot(4,1,2); imagesc( tempCrossings' );
        sp(3) = subplot(4,1,3); imagesc( tempEvents' )
        linkaxes(sp, 'xy');
        impixelinfo;       
    end

    for roi = 1:expt.Nroi
        connComp = bwconncomp( tempEvents(:,roi) );
        Nevent(run,roi) = connComp.NumObjects;
        for evt = 1:Nevent(run,roi)
            fluorEvent{run,roi}(evt).periScan = connComp.PixelIdxList{evt};
            localEvtScan = find(fluor(run).z.ROI(fluorEvent{run,roi}(evt).periScan, roi) > zThresh);
            fluorEvent{run,roi}(evt).eventScan = fluorEvent{run,roi}(evt).periScan(localEvtScan);
            fluorEvent{run,roi}(evt).NeventScan = numel(fluorEvent{run,roi}(evt).eventScan);
            fluorEvent{run,roi}(evt).Tperi = T{run}(fluorEvent{run,roi}(evt).periScan);
            fluorEvent{run,roi}(evt).z = fluor(run).z.ROI(fluorEvent{run,roi}(evt).periScan, roi);
            fluorEvent{run,roi}(evt).zMax = max( fluorEvent{run,roi}(evt).z );
            fluorEvent{run,roi}(evt).dFF = fluor(run).dFF.ROI(fluorEvent{run,roi}(evt).periScan, roi);
            fluorEvent{run,roi}(evt).dFFmax = max( fluorEvent{run,roi}(evt).dFF );
            fluorEvent{run,roi}(evt).dFFbase = prctile( fluorEvent{run,roi}(evt).dFF, 10); %max( fluor(r).dFF.ROI(fluorEvent{r,roi}(evt).periScan, roi) );
            %{
            if show
                subplot(4,1,4);
                plot( fluorEvent{run,roi}(evt).Tperi, fluorEvent{run,roi}(evt).dFF, 'k' ); hold on;
                plot( fluorEvent{run,roi}(evt).Tperi(localEvtScan), fluorEvent{run,roi}(evt).dFF(localEvtScan), 'bo' ); hold on;
                xlabel('Peri-Event Time (s)'); ylabel('Fluor'); 
                title( sprintf('%s [run, ROI, event] = [%i, %i, %i] ', expt.name, run, roi, evt ), 'Interpreter','none' );
                pause;
                cla;
            %}
        end
    %}