% Find fluorEvent associated with stillness and locomotion - run BoutResponse first
windowPad = [-2,2]; padDur = diff(windowPad);
minIso = 10;
SAFE = cell(1,Nexpt); Nsafe = cell(1,Nexpt); LAFE = cell(1,Nexpt); Nlafe = cell(1,Nexpt); CAFE = cell(1,Nexpt); Ncafe = cell(1,Nexpt);
tic; close all;
for x = xPresent
    [stillEpoch{x}, stillSumm{x}] = GetStillEpochs(expt{x}, Tscan{x}, loco{x}, fluor{x}, deform{x}, csdBout{x}, 'show',true, 'criterion','state'); % find long periods of stillness between locomotive bouts

    SAFE{x} = cell(1,expt{x}.Nroi);  LAFE{x} = cell(1,expt{x}.Nroi);
    stillSumm{x}.Nevent = nan(stillSumm{x}.Nepoch, expt{x}.Nroi);  stillSumm{x}.eventRate = nan(stillSumm{x}.Nepoch, expt{x}.Nroi);
    stillSumm{x}.totEventDur = nan(stillSumm{x}.Nepoch, expt{x}.Nroi); stillSumm{x}.durRatio = nan(stillSumm{x}.Nepoch, expt{x}.Nroi);
    LAFEsumm(x).boutDur = [];
    LAFEsumm(x).magPeak = cell(1,expt{x}.Nroi); LAFEsumm(x).totEventDur = cell(1,expt{x}.Nroi); LAFEsumm(x).totBoutDur = cell(1,expt{x}.Nroi);  LAFEsumm(x).durRatio = cell(1,expt{x}.Nroi); 
    LAFEsumm(x).startDelay = cell(1,expt{x}.Nroi); LAFEsumm(x).peakDelay = cell(1,expt{x}.Nroi); 
    for run = expt{x}.preRuns
        periBout{x}(run).Nevent = nan(periBout{x}(run).Nbout, expt{x}.Nroi);
        periBout{x}(run).eventRate = nan(periBout{x}(run).Nbout, expt{x}.Nroi);
        periBout{x}(run).bIso = find(periBout{x}(run).iso(:,1) > minIso)';
        periBout{x}(run).raster = cell(1,periBout{x}(run).Nbout);  %#ok<*SAGROW>
        for b = periBout{x}(run).bIso
            periBout{x}(run).raster{b} = false( size(periBout{x}(run).fluor{b}) );
        end
        LAFEsumm(x).boutDur = [LAFEsumm(x).boutDur, periBout{x}(run).dur(periBout{x}(run).bIso)+padDur]; % duration (including pad time) of each isolated bout
        for roi = find(cellfun(@numel, fluorEvent{x}{run} ) > 0) %1:expt{x}.Nroi % 1:expt{x}.Nruns
            tempTstart = [fluorEvent{x}{run}{roi}.Tstart];
            % Find fluor events within epochs of stillness
            for s = find([stillEpoch{x}.run] == run) %stillSumm{x}.sPre %1:stillSumm{x}.Nepoch
                eStill = find(tempTstart >= stillEpoch{x}(s).Tstart & tempTstart <= stillEpoch{x}(s).Tstop);
                if ~isempty(eStill)
                    tempSAFE = fluorEvent{x}{run}{roi}(eStill);
                    SAFE{x}{roi} = [SAFE{x}{roi}, tempSAFE];
                    stillSumm{x}.Nevent(s,roi) = numel(eStill); % stillEpoch{x}(s).Nevent(roi)
                    stillSumm{x}.eventRate(s,roi) = 60*stillSumm{x}.Nevent(s,roi)/stillEpoch{x}(s).roiDur(roi); % events per minute %stillEpoch{x}(s).eventRate(roi)
                    stillSumm{x}.totEventDur(s,roi) = sum([tempSAFE.dur]);
                    stillSumm{x}.durRatio(s,roi) = stillSumm{x}.totEventDur(s,roi)/stillEpoch{x}(s).dur; % fraction of time spent active  sum([tempSAFE.dur])
                    for e = eStill
                        [~,~,tempScan] = intersect( fluorEvent{x}{run}{roi}(e).T, stillEpoch{x}(s).T );
                        stillEpoch{x}(s).raster(tempScan,roi) = true;  
                    end
                else
                    stillSumm{x}.Nevent(s,roi) = 0;
                    stillSumm{x}.eventRate(s,roi) = 0;
                    stillSumm{x}.totEventDur(s,roi) = 0;
                    stillSumm{x}.durRatio(s,roi) = 0;
                end
            end
            % Find fluor events around locomotive bouts
            for b = periBout{x}(run).bIso %1:periBout{x}(run).Nbout
                periBoutWindow = [periBout{x}(run).T{b}(periBout{x}(run).boutScan{b}(1)),  periBout{x}(run).T{b}(periBout{x}(run).boutScan{b}(end))] + windowPad;
                eLoco = find(tempTstart >= periBoutWindow(1) & tempTstart <= periBoutWindow(2));
                LAFE{x}{roi} = [LAFE{x}{roi}, fluorEvent{x}{run}{roi}(eLoco)];
                periBout{x}(run).Nevent(b,roi) = numel(eLoco);
                periBout{x}(run).eventRate(b,roi) = 60*periBout{x}(run).Nevent(b,roi)/(periBout{x}(run).dur(b)+padDur); % events per minute
                if periBout{x}(run).Nevent(b,roi) > 0
                    for e = eLoco
                        [~,~,tempScan] = intersect( fluorEvent{x}{run}{roi}(e).T, periBout{x}(run).T{b} );
                        periBout{x}(run).raster{b}(tempScan,roi) = true;  % periBout{x}(run).event{b}(E,roi).boutScan
                    end
                    tempLAFE = fluorEvent{x}{run}{roi}(eLoco); %[periBout{x}(run).event{b}(1:E,roi)];
                    %LAFEsumm(x).eventTime = sum([tempLAFE.dur]);
                    startDelay = tempLAFE(1).Tstart - periBout{x}(run).Tstart(b);
                    LAFEsumm(x).startDelay{roi} = [LAFEsumm(x).startDelay{roi}, startDelay];
                    tempTotDur = sum([tempLAFE.dur]);
                    LAFEsumm(x).totEventDur{roi} = [LAFEsumm(x).totEventDur{roi}, tempTotDur];
                    LAFEsumm(x).durRatio{roi} = [LAFEsumm(x).durRatio{roi}, tempTotDur/(periBout{x}(run).dur(b)+padDur)];
                    [maxPeak, ePeak] = max([tempLAFE.magPeak]);
                    LAFEsumm(x).magPeak{roi} = [LAFEsumm(x).magPeak{roi}, maxPeak];
                    peakDelay = tempLAFE(ePeak).Tpeak - periBout{x}(run).Tstart(b);
                    LAFEsumm(x).peakDelay{roi} = [LAFEsumm(x).peakDelay{roi}, peakDelay]; % tempLAFE(ePeak).peakDelay
                else
                    LAFEsumm(x).totEventDur{roi} = [LAFEsumm(x).totEventDur{roi}, 0];
                    LAFEsumm(x).durRatio{roi} = [LAFEsumm(x).durRatio{roi}, 0];
                end
            end
        end
    end
    if stillSumm{x}.Nepoch > 0 
        stillSumm{x}.totEventRate = 60*sum([stillSumm{x}.Nevent], 1, 'omitnan')./sum(stillSumm{x}.roiDur, 1); % total events per minute total ROI time
    else
        stillSumm{x}.totEventRate = nan(1, expt{x}.Nroi);
    end
    %{
    if ~isnan(expt{x}.csd)
        CAFE{x} = cell(1,expt{x}.Nroi);
        for run = expt{x}.csd           
            csdWindow = [csdBout{x}(run).Tstart, csdBout{x}(run).Tstop] + windowPad;
            for roi = find(NfluorEvent{x}(run,:)) %1:expt{x}.Nroi
                tempTstart = [fluorEvent{x}{run}{roi}.Tstart];
                eCSD = find(tempTstart >= csdWindow(1) & tempTstart <= csdWindow(2));
                CAFE{x}{roi} = [CAFE{x}{roi}, fluorEvent{x}{run}{roi}(eCSD)];
            end
        end
        Ncafe{x} = cellfun(@numel, CAFE{x});
    end
    %}
    Nsafe{x} = cellfun(@numel, SAFE{x});  Nlafe{x} = cellfun(@numel, LAFE{x}); 
    % Pool runs and summarize LAFE results
    LAFEsumm(x).Nevent = vertcat(periBout{x}.Nevent); % # of fluor events for each qualified bout (row) and ROI (column)
    LAFEsumm(x).Nevent(isnan(LAFEsumm(x).Nevent(:,1)),:) = [];
    LAFEsumm(x).Nbout = size(LAFEsumm(x).Nevent,1); % # of qualified bouts
    LAFEsumm(x).eventRate = vertcat( periBout{x}(:).eventRate );
    LAFEsumm(x).eventRate(isnan(LAFEsumm(x).eventRate(:,1)),:) = [];
    %LAFEsumm(x).rateDiff = median(LAFEsumm(x).eventRate,1,'omitnan') - median(stillSumm{x}.eventRate,1,'omitnan');
    LAFEsumm(x).ratioDiff = cellfun(@median, LAFEsumm(x).durRatio) - median(stillSumm{x}.durRatio, 1, 'omitnan'); 
    LAFEsumm(x).responder = []; LAFEsumm(x).responder_pos = []; LAFEsumm(x).responder_neg = [];
    for roi = 1:expt{x}.Nroi
        try
            pROIrate = ranksum(stillSumm{x}.durRatio(:,roi), LAFEsumm(x).durRatio{roi}' ); %JitterPlot( cell2padmat({stillSumm{x}.durRatio(:,roi), LAFEsumm(x).durRatio{roi}}) )
            if pROIrate < 0.05 && LAFEsumm(x).ratioDiff(roi) ~= 0 %  LAFEsumm(x).rateDiff(roi)
                LAFEsumm(x).responder = [LAFEsumm(x).responder, roi]; 
                %JitterPlot( cell2padmat({stillSumm{x}.eventRate(:,roi), LAFEsumm(x).eventRate(:,roi)}) )
            end
        catch
            fprintf('\nranksum failed for [x,roi] = [%i, %i]', x, roi);
        end
    end
    LAFEsumm(x).responder_pos =  LAFEsumm(x).responder( LAFEsumm(x).ratioDiff(LAFEsumm(x).responder) > 0 );    
    LAFEsumm(x).responder_neg =  LAFEsumm(x).responder( LAFEsumm(x).ratioDiff(LAFEsumm(x).responder) < 0 );
    LAFEsumm(x).Nresponder = numel(LAFEsumm(x).responder);
    LAFEsumm(x).Npos = numel(LAFEsumm(x).responder_pos);
    LAFEsumm(x).Nneg = numel(LAFEsumm(x).responder_neg);   
    LAFEsumm(x).boutRespFrac = sum(LAFEsumm(x).Nevent>0, 2)/expt{x}.Nroi; %LAFEsumm(x).Nbout; % Fraction of units that responded to each bout
    LAFEsumm(x).boutActFrac = sum(LAFEsumm(x).boutRespFrac>0, 1)/LAFEsumm(x).Nbout; % Fraction of bouts that activated at least one unit
    LAFEsumm(x).roiActFrac = sum(LAFEsumm(x).Nevent>0, 1)/LAFEsumm(x).Nbout; % fraction of bouts that activated each unit
    %LAFEsumm(x).responder = find(LAFEsumm(x).roiActFrac > 0.5); % identify ROI that responded to at least 50% of bouts
    LAFEsumm(x).responderFrac = LAFEsumm(x).Nresponder/expt{x}.Nroi; % fraction of ROI that responded to at least 50% of bouts
    LAFEsumm(x).posFrac = LAFEsumm(x).Npos/expt{x}.Nroi; 
    LAFEsumm(x).negFrac = LAFEsumm(x).Nneg/expt{x}.Nroi; 
    LAFEsumm(x).startDelayMed = cellfun(@median, LAFEsumm(x).startDelay);
    LAFEsumm(x).startDelayMed(setdiff(1:expt{x}.Nroi, LAFEsumm(x).responder_pos)) = NaN; % exclude non-positive responders
    LAFEsumm(x).peakDelayMed = cellfun(@median, LAFEsumm(x).peakDelay);
    LAFEsumm(x).peakDelayMed(setdiff(1:expt{x}.Nroi, LAFEsumm(x).responder_pos)) = NaN;
    LAFEsumm(x).totDurMed = cellfun(@median, LAFEsumm(x).totEventDur);
    LAFEsumm(x).totDurMed(setdiff(1:expt{x}.Nroi, LAFEsumm(x).responder_pos)) = NaN;
    LAFEsumm(x).durRatioMed = cellfun(@median, LAFEsumm(x).durRatio);
    LAFEsumm(x).durRatioMed(setdiff(1:expt{x}.Nroi, LAFEsumm(x).responder_pos)) = NaN;
    LAFEsumm(x).magPeakMed = cellfun(@median, LAFEsumm(x).magPeak);
    LAFEsumm(x).magPeakMed(setdiff(1:expt{x}.Nroi, LAFEsumm(x).responder_pos)) = NaN;
    toc
    
    % Cluster ROI into putative axons
    %[axon{x}, expt{x}, stillCosSim{x}, stillCorr{x}] = MergeROI3D(expt{x}, Tscan{x}, loco{x}, ROI{x}, fluor{x}, SAFE{x}, 'method','cluster', 'mergeThresh',1.5, 'show',false); 
end

%% Show epoch data vs raster
close all;
figure('Units','normalized', 'OuterPosition', [0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
for x = xPresent
    for s = stillSumm{x}.sPre
        Ttick = 0:5:stillEpoch{x}(s).dur;
        Xtick = round(expt{x}.scanRate*Ttick);
        
        
        sp(1) = subplot(2,1,1);
        imagesc( stillEpoch{x}(s).dFF' )
        set(gca,'Xtick',Xtick', 'TickDir','out');
        
        sp(2) = subplot(2,1,2);
        imagesc( stillEpoch{x}(s).raster' );
        

        set(gca,'Xtick',Xtick', 'XtickLabel',Ttick, 'TickDir','out');
        xlabel('Epoch Time (s)');
        linkaxes(sp,'xy');
        impixelinfo;
        pause;
        clf;
    end
end

%% Show peribout data vs raster
tickTime = 2; % interval, in seconds, to use for x axis ticks
close all;
figure('Units','normalized', 'OuterPosition', [0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
for x = xPresent
    tickInt = round(tickTime*expt{x}.scanRate);
    for run = expt{x}.preRuns
        for b = periBout{x}(run).bIso
            Tbout = periBout{x}(run).T{b} - periBout{x}(run).Tstart(b);
            zeroTick = find(Tbout==0);
            Xtick = unique([flip(zeroTick:-tickInt:1), zeroTick:tickInt:numel(Tbout)]);
            Ttick = round(Tbout(Xtick));
            sp(1) = subplot(2,1,1);
            imagesc( periBout{x}(run).fluor{b}' )
            set(gca,'Xtick',Xtick, 'XtickLabel',Ttick, 'TickDir','out');
            title( sprintf('[x,run,b] = [%i, %i, %i]', x, run, b) );

            sp(2) = subplot(2,1,2);
            imagesc( periBout{x}(run).raster{b}' );
            set(gca,'Xtick',Xtick, 'XtickLabel',Ttick, 'TickDir','out');
            xlabel('Peri-bout Time (s)');
            linkaxes(sp,'xy');
            impixelinfo;
            pause;
            clf;
        end
    end
end

%%  Characterize spontaneous activity during quiet wakefulness, broadly

% x% of ROI exhibited activity during quiet wakefulness
% median duration of event
% median event rate
tempStillSumm = [stillSumm{xPresent}];
totEventRatePool = [tempStillSumm.totEventRate];
ecdf(totEventRatePool)
median(totEventRatePool)
SAFEpool = [];
for x = xPresent
   SAFEpool = [SAFE{x}{:}];
end
median([SAFEpool.dur])