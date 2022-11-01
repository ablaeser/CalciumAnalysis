function fluor = GetAxonFluor(expt, sbxInfo, axon, fluor, deform, varargin) 
%
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'axon', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'deform', @isstruct )
addParameter( IP, 'window', 101, @isnumeric) %rolling percentile window
addParameter( IP, 'basePrct', 10, @isnumeric) %rolling percentile value
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, sbxInfo, axon, fluor, deform, varargin{:} );  %
windowSize = IP.Results.window;
if mod(windowSize,2) ==0,  windowSize = windowSize+1; end
basePrct = IP.Results.basePrct; 
overwrite = IP.Results.overwrite;
Naxon = numel(axon);
tic
if ~isfield(fluor(1).F, 'axon') || overwrite
    savePath = strcat( expt.dir, expt.name, '_axonFluor.mat' );
    if ~exist(savePath, 'file') || overwrite
        Faxon = nan(expt.totScan, Naxon); Fnp = nan(expt.totScan, Naxon);
        fprintf('\nExtracting traces');      
        w = waitbar(0,'Extracting axonal and neuropil traces...'); % w = 
        for s = 1:expt.totScan-1
            tempVol = double(readSBX(expt.sbx.interp, sbxInfo, s, 1, -1, [])); % %double( readSBX(expt.sbx, sbxInfo, expt.Nplane*(s-1)+1, expt.Nplane, 1, []) ); 
            for a = 1:Naxon
                Faxon(s,a) = mean( tempVol(axon(a).ind) );
                Fnp(s,a) = mean( tempVol(axon(a).neuropil) );
            end
            waitbar(s/(expt.totScan-1));
        end
        clearvars tempVol; delete(w); 
        F = Faxon - Fnp + mean(Fnp,1,'omitnan'); % plot( F(:,1) );
        toc

        % Censor scans where each ROI's planes had a bad deformation value
        transAPcat = cat(1, deform.transAP);
        for a = 1:Naxon
            nan_scan = find( isnan(sum(transAPcat(:,axon(a).zUse),2) ) )';
            F(nan_scan,a) = NaN;
            fprintf('\nROI %i: censored %i scans with bad deformation values', a, numel(nan_scan) );
        end

        % Use the censored, neuropil-subtracted signal to normalize
        Fo = MovingPercentile(F, basePrct, windowSize, 'pre'); % movprctile(F, basePrct, windowSize, 1); %plot([F(:,1), Fbase(:,1)] );
        dFF = (F-Fo)./Fo;

        % Deconvolve dF/F
        activity = nan(size(dFF,1),size(dFF,2));
        %{
        if expt.Nplane == 1
            w = waitbar(0, 'Deconvolving censored dF/F');
            for a = flip(1:Naxon)
                try
                    goodFrames = find(~isnan(dFF(:,a)))';
                    dFFcens = dFF(goodFrames,a);
                    [~, actCens, ~] = deconvolveCa( dFFcens, 'exp2', 'foopsi', 'optimize_pars', 'optimize_smin', 'maxIter',20, 'window',300 );  %ar1  deconvOpts(r)
                    activity(goodFrames,a) = actCens;
                    %figure; yyaxis left; plot( dFF(:,r) ); yyaxis right; plot( activity(:,r) ); 
                catch
                    fprintf('\nROI %i:  Deconvolution failed', a);
                end
                waitbar(a/Naxon, w);
            end
            delete(w);
        end
        %}
        toc
        fprintf('\nSaving %s', savePath);
        save(savePath, 'Faxon','Fnp','F','Fo','dFF','activity','nan_scan', '-v7.3');
    else
        fprintf('\nLoading %s', savePath);
        load(savePath, 'Faxon','Fnp','F','Fo','dFF','activity' ); % ,'nan_scan'
    end

    % Break fluor signals up by run
    fprintf('\nBreaking fluor signals up by run...  ')
    scanLims = expt.scanLims;
    for runs = 1:expt.Nruns
        fluor(runs).Froi.axon = Faxon(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Fnp.axon = Fnp(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).F.axon = F(scanLims(runs)+1:scanLims(runs+1),:); 
        fluor(runs).Fo.axon = Fo(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).dFF.axon = dFF(scanLims(runs)+1:scanLims(runs+1),:);
        fluor(runs).z.axon = normalize(fluor(runs).dFF.axon, 1); % 
        fluor(runs).act.axon = activity(scanLims(runs)+1:scanLims(runs+1),:);
        % Correlation between axons for this run
        fluor(runs).Froi.corr = corr( fluor(runs).Froi.axon, 'Rows','complete' );
        fluor(runs).Fnp.corr = corr( fluor(runs).Fnp.axon, 'Rows','complete' );
        fluor(runs).F.corr = corr( fluor(runs).F.axon, 'Rows','complete' );
        fluor(runs).Fo.corr = corr( fluor(runs).Fo.axon, 'Rows','complete' );
        fluor(runs).dFF.corr = corr( fluor(runs).dFF.axon, 'Rows','complete' );
        fluor(runs).act.corr = corr( fluor(runs).act.axon, 'Rows','complete' );
        %{
        for a = 1:Naxon
            subplot(4,2,[1,3,5,7]);
            %imshow( axon(a).labelFoot, [] );
            imshow( label2rgb(axon(a).labelFoot), [] );
            hold on;
            title(sprintf('[run, axon] = [%i, %i]', runs, a));
            impixelinfo
            
            sp(1) = subplot(4,2,2); cla;
            plot( fluor(runs).F.ROI(:,axon(a).ROI) ); % hold on;
            ylabel('F'); title('ROI');
            
            sp(2) = subplot(4,2,4);
            plot( fluor(runs).F.axon(:,a) ); 
            ylabel('F'); title('Axon');            
            
            
            sp(3) = subplot(4,2,6); cla;
            plot( fluor(runs).dFF.axon(:,a) ); ylabel('dF/F');
            
            sp(4) = subplot(4,2,8); cla;
            plot( fluor(runs).act.axon(:,a) ); ylabel('Activity');
            
            linkaxes(sp,'x');
            xlim([-Inf,Inf]);
            pause;
        end
        %}
    end
    toc
else
    fprintf('\nROI fluor already found!');
end