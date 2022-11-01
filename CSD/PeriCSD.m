function [csdBout, csdStat, csdParam] = PeriCSD(expt, T, loco, deform, fluor, defVars, varargin) % , periParam
IP = inputParser;
addRequired( IP, 'expt', @isstruct ) %addRequired( IP, 'infoStruct', @isstruct )
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'loco', @isstruct )
addRequired( IP, 'deform', @isstruct )
addRequired( IP, 'fluor', @isstruct )
addRequired( IP, 'defVars', @iscell )
addParameter( IP, 'window', 300, @isnumeric ); % peak must occur within this many seconds of starting movie
addParameter( IP, 'thresh', 0.1, @isnumeric ) % fraction of peak used to determine limits of CSD
addParameter( IP, 'iso', 5, @isnumeric )
addParameter( IP, 'base', 10, @isnumeric )
addParameter( IP, 'run', 4, @isnumeric )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, T, loco, deform, fluor, defVars, varargin{:} );  % , periParam
csdWindowTime = IP.Results.window;
threshFrac = IP.Results.thresh;
csdParam.base = IP.Results.base;
csdParam.run = IP.Results.run;
csdParam.iso = IP.Results.iso;
csdParam.min_vel_on = -Inf;
overwrite = IP.Results.overwrite;
show = IP.Results.show;

savePath = strcat( expt.dir, expt.name, '_CSDbout.mat' );
if exist(savePath, 'file') && ~overwrite
    fprintf('\nLoading %s', savePath);
    load(savePath, 'csdBout', 'csdStat', 'csdParam', 'threshInd', 'cutInd', 'peakInd', 'Ffilt', 'dFfilt', 'FcsdPre');
else
    NdefVars = numel(defVars);
    dT = 1/expt.scanRate;
    if isempty(csdParam.iso), csdParam.iso = csdParam.base; end
    csdParam.Nshuff = 1;
    csdParam.sigThresh = 0.05;
    csdParam.on = 0; %csdParam.run;
    csdParam.dT = 1/expt.scanRate;
    csdParam.NbaseScan = ceil(expt.scanRate*csdParam.base);
    csdParam.NrunScan = ceil(expt.scanRate*csdParam.run);
    csdParam.T = csdParam.dT*(-csdParam.NbaseScan:csdParam.NrunScan);
    
    % Use derivative for an initial estimate of the edges of CSD wave
    gaussWidth = 5; gaussSigma = 3;
    gaussFilt = MakeGaussFilt( gaussWidth, 0, gaussSigma, expt.scanRate );
    Ffilt = filtfilt(gaussFilt,1,fluor.F.vol);
    dFfilt = diff( Ffilt );
    Nscan = size(Ffilt,1);
    csdWindowScans = 1:min(ceil(csdWindowTime/dT), Nscan-1 ); % expt.Nscan
    [~,cutInd(1)] = max( dFfilt(csdWindowScans) );
    [~,cutInd(2)] = min( dFfilt(cutInd(1):csdWindowScans(end)) );
    cutInd(1) = cutInd(1) + 1;
    cutInd(2) = cutInd(2) + cutInd(1) + 1; % derivative drops one index
    
    % Use peak to define CSD
    FcsdPre = prctile( fluor.F.vol(1:cutInd(1)-1), 10 );
    dFFcsd = (fluor.F.vol-FcsdPre)./FcsdPre;
    [dFFcsdPeak, peakInd] = max( dFFcsd(cutInd(1):cutInd(2)) );
    peakInd = peakInd + cutInd(1);
    dFFthresh = threshFrac*dFFcsdPeak;
    threshInd(1) = find( dFFcsd(1:peakInd-1) < dFFthresh, 1, 'last' ); %
    postThreshInd = find( dFFcsd(peakInd:end) < dFFthresh, 1, 'first' );
    if isempty(postThreshInd)
        threshInd(2) = cutInd(2);
    else
        threshInd(2) = postThreshInd + peakInd;
    end
    
    % Grab peri-CSD data and summarize
    csdScans = threshInd(1):threshInd(2);
    csdBout = GetBoutData( {csdScans}, T, loco, deform, fluor.F.ROI, defVars, csdParam ); % fluor.z.ROI
    csdStat = struct();  %'speed',[], 'trans_x',[], 'trans_y',[], 'scale_x',[], 'scale_y',[], 'shear_x',[], 'shear_y',[], 'shift_z',[], 'dTransR',[], 'stretch_x',[], 'stretch_y',[], 'dShiftZ',[], 'fluor',[]
    if csdBout.Nbout > 0
        % Summarize effect on readout variables before/during/after bout
        csdStat.speed = BoutEffect( csdBout.speed, csdBout.preScan, csdBout.boutScan, csdBout.postScan );
        csdStat.fluor = BoutEffect( csdBout.fluor, csdBout.preScan, csdBout.boutScan, csdBout.postScan );
        for v = 1:NdefVars
            csdStat.(defVars{v}) = BoutEffect( csdBout.(defVars{v}), csdBout.preScan, csdBout.boutScan, csdBout.postScan );
        end
        % Compare effect sizes to those obtained from shuffled data
        csdStat.speed.pEffect = BoutEffectPval(csdBout, csdStat, csdParam, 'speed');
        csdStat.fluor.pEffect = BoutEffectPval(csdBout, csdStat, csdParam, 'fluor');
        for v = 1:NdefVars
            csdStat.(defVars{v}).pEffect = BoutEffectPval(csdBout, csdStat, csdParam, defVars{v});
        end
    end
    
    % Save the results
    fprintf('\nSaving %s', savePath);
    save(savePath, 'csdBout', 'csdStat', 'csdParam', 'threshInd', 'cutInd', 'peakInd', 'Ffilt', 'dFfilt', 'FcsdPre', '-v7.3'); 
end

if show
    Tcsd = T-T(threshInd(1));
    %close all;
    figure('WindowState','maximized', 'Color','w');
    sp(2) = subplot(2,1,2);
    plot( Tcsd(2:end), dFfilt ); hold on;
    plot( Tcsd(cutInd), dFfilt(cutInd-1), 'ro' );
    set(gca, 'TickDir','out', 'tickLength',[0.01,0], 'box','off')
    ylabel('Time Derivative (filtered)'); xlabel('Peri-CSD Time (s)');
    
    sp(1) = subplot(2,1,1);
    plot( Tcsd, fluor.F.vol ); hold on;
    plot( Tcsd, Ffilt );
    plot( Tcsd(cutInd), fluor.F.vol(cutInd), 'ro' )
    plot( Tcsd(threshInd), fluor.F.vol(threshInd), 'ko' )
    plot( Tcsd(peakInd), fluor.F.vol(peakInd), 'kx');
    line( Tcsd(threshInd), FcsdPre*[1,1], 'lineStyle','-', 'color','k' )
    set(gca, 'TickDir','out', 'tickLength',[0.01,0], 'box','off')
    ylabel('dF/Fo'); % Fluor
    title( sprintf('%s', expt.name ), 'Interpreter','none' );
    xlim( [-Inf,Inf] );
    linkaxes(sp,'x');
    %{
    plot( Tcsd, dFFcsd ); hold on;
    plot( Tcsd(peakInd), dFFcsd(peakInd), 'o' );
    line( Tcsd([1, numel(dFFcsd)]), dFFcut*[1,1], 'lineStyle','--', 'color','r');
    line( Tcsd(cutInd), [0,0], 'lineStyle','-', 'color','k' )
    set(gca, 'TickDir','out', 'tickLength',[0.01,0], 'box','off')
    ylabel('dF/Fo'); % Fluor
    xlabel('Time (s)');
    xlim( [-Inf,Inf] ); % xlim( Tcsd(csdInd([1,end])) );
    %}
end
end