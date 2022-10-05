function [maxLikeTheta, likeTheta, errStruct, P, Pdens, trainData, testData, decodeLim] = DecodeDirection(deformIn, activityIn, varargin)
% Decode the direction of stretch from the neural activity using the maximum likelihood method
IP = inputParser;
addRequired( IP, 'stretchIn', @isnumeric )
addRequired( IP, 'activityIn', @isnumeric )
addParameter( IP, 'minFrac', 0.5, @isnumeric )
addParameter( IP, 'trainFrac', 75, @isnumeric )
addParameter( IP, 'thetaLim', [-pi,pi], @isnumeric )
addParameter( IP, 'dTheta', pi, @isnumeric )
addParameter( IP, 'thetaRef', 0, @isnumeric )
addParameter( IP, 'thetaShift', 0, @isnumeric )
addParameter( IP, 'NactBin', 30, @isnumeric )
addParameter( IP, 'Npseudo', 1, @isnumeric )
addParameter( IP, 'Nboot', 1000, @isnumeric )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'pause', Inf, @isnumeric )
parse( IP, deformIn, activityIn, varargin{:} ); 
minFrac = IP.Results.minFrac;
trainFrac = IP.Results.trainFrac;
thetaLim = IP.Results.thetaLim;
dTheta = IP.Results.dTheta;
thetaRef = IP.Results.thetaRef;
thetaShift = IP.Results.thetaShift;
NactBin = IP.Results.NactBin;
Npseudo = IP.Results.Npseudo;
Nboot = IP.Results.Nboot;
show = IP.Results.show;
subName = IP.Results.name;
pauseTime = IP.Results.pause;

tic;
Nroi = size(activityIn, 2);
% Restrict data used for decoding by stretch angle and magnitude
fprintf('\n  Restricting data... ');
NinFrames = size(deformIn, 1);
[ dirIn, ~, magIn, binIn, ~, ~, statsIn ] = AngularDist( deformIn(:,1), deformIn(:,2), 'thetaRef',thetaRef, 'dThetaDir',dTheta, 'shift',thetaShift, 'dPct',5, 'show',false); %true
goodAngleBins = find( binIn.dir.cent >= thetaLim(1) & binIn.dir.cent <= thetaLim(2) )';
decodeLim = struct('dir',thetaLim, 'mag',nan(1,2), 'frames',[]);
decodeLim.mag(2) = min( statsIn.dir.pct(goodAngleBins,end-1) );
decodeLim.mag(1) = minFrac*decodeLim.mag(2); % decodeLim.mag(2)/2;
decodeLim.frames = find( dirIn >= thetaLim(1) & dirIn <= thetaLim(2) & magIn >= decodeLim.mag(1) & magIn <= decodeLim.mag(2) );
decodeStretch = deformIn(decodeLim.frames,:);
activityIn(activityIn<0) = 0; % negative activity values cause problems using log
decodeActivity = activityIn(decodeLim.frames,:);
%decodeDir = decodeDir(decodeLim.frames,:);
%decodeMag = decodeMag(decodeLim.frames,:);
fprintf('found %i good frames', numel(decodeLim.frames) );

% Define training and testing data subsets
fprintf('\n  Defining training/testing subsets: ');
trainData = struct('N',NaN, 'frames',[], 'dir',[], 'mag',[], 'bin',[] );
testData = struct('N',NaN, 'frames',[], 'dir',[], 'mag',[], 'bin',[], 'trueBin',[] );
% Bin the restricted data set by angle, then divide into training/testing subsets, with equal numbers of training samples for all bins
angBinData = struct('N',NaN, 'frames',[], 'dir',[], 'mag',[], 'bin',[] );
[angBinData.dir, ~, angBinData.mag, angBinData.bin, ~, ~, ~] = AngularDist( decodeStretch(:,1), decodeStretch(:,2), 'thetaRef',thetaRef, 'dThetaDir', dTheta, 'show',false ); % use the merged versions of stretch, to be consistent with GLM. By convention consider ML to be the x-axis
binN = cellfun( @numel, angBinData.bin.dir.ind );
trainN = round((trainFrac/100)*min( binN( binN > 0 ) ));
trainData.bin = angBinData.bin;
trueTestBinCell = cell(1, angBinData.bin.dir.Nbin );
bUse = find( binN > 0 ); %find(~cellfun( @isempty, trainData.bin.dir.ind ) ); 
Nuse = numel(bUse);
for b = bUse
    binTrainInd = randsample( binN(b), trainN );
    binTestInd = 1:binN(b); 
    binTestInd( binTrainInd ) = [];
    trainData.bin.dir.ind{b} = angBinData.bin.dir.ind{b}( binTrainInd );
    testData.bin.dir.ind{b} = angBinData.bin.dir.ind{b}( binTestInd );
    trueTestBinCell{b} = repmat( b, numel( testData.bin.dir.ind{b} ), 1 );
end
trainData.frames = vertcat( trainData.bin.dir.ind{:} );
trainData.N = numel( trainData.frames );
trainData.mag = angBinData.mag( trainData.frames, :);
trainData.dir = angBinData.dir( trainData.frames, :);
trainData.activity = decodeActivity(trainData.frames,:); 
testData.frames = vertcat( testData.bin.dir.ind{:} );
testData.N = numel( testData.frames );
testData.dir = angBinData.dir(testData.frames,:);
[testData.dir, testDirSortInd] = sort( angBinData.dir(testData.frames,:), 'ascend');
testData.frames = testData.frames(testDirSortInd); % sort frames by direction
testData.mag = angBinData.mag(testData.frames,:);
testData.activity = decodeActivity(testData.frames,:); 
testData.bin.dir.Nbin = angBinData.bin.dir.Nbin;
testData.bin.dir.lim = angBinData.bin.dir.lim;
testData.bin.dir.cent = angBinData.bin.dir.cent;
testData.trueBin = vertcat( trueTestBinCell{:} );

fprintf('found %i training frames (%i per bin),  %i testing frames', trainData.N, trainData.N/Nuse, testData.N );
%testData.bin.dir.ind(bUse)

% Show the training/testing data
%{
if show
    figure('color','w');
    polarplot( dirIn, magIn, 'k.' ); hold on; 
    polarplot( trainData.dir, trainData.mag, 'b.' ); 
    polarplot( testData.dir, testData.mag, 'r.' );
    title('Deform Samples'); 
    legend('All','Train','Test', 'Location','best');
end
%}
% Determine activity probability densities from training data
fprintf('\n  Training the model...')
Pdens = struct('a',[], 'a_theta',[], 'binEdges',[]);
[~,Pdens.binEdges] = histcounts( log10(trainData.activity(:)), NactBin ); 
Pdens.binEdges = [-Inf, Pdens.binEdges, Inf];
Pdens.a = cell(1,Nroi); Pdens.a_theta = cell(trainData.bin.dir.Nbin, Nroi);   
for roi = 1:Nroi
    tempCounts = histcounts( log10(trainData.activity(:,roi)), Pdens.binEdges, 'Normalization','count'); % fitdist( trainData.activity(:,r), 'Exponential'  );
    tempCounts = tempCounts + Npseudo; % add pseudocount to prevent P == 0
    Pdens.a{roi} = tempCounts/sum(tempCounts);       
    for b = flip(1:trainData.bin.dir.Nbin)
        tempCounts = histcounts( log10(decodeActivity(trainData.bin.dir.ind{b},roi)), Pdens.binEdges, 'Normalization','count'); % new version
        %tempCounts = histcounts( log10(trainData.activity(trainData.bin.dir.ind{b},r)), Pdens.binEdges, 'Normalization','count'); %old version
        tempCounts = tempCounts + Npseudo; % prevent P = 0
        Pdens.a_theta{b,roi} = tempCounts/sum(tempCounts); %fitdist( trainData.activity(trainData.bin.dir.ind{b},r), 'Exponential'  );
    end    
end

% Get probability of theta from each test frame's activity:  
fprintf('\n  Calculating probabilities...')
P = struct('theta',nan(1, 1, trainData.bin.dir.Nbin), 'a',nan(testData.N, Nroi), 'a_theta',nan(testData.N, Nroi, testData.bin.dir.Nbin ), 'theta_a',[] );
for b = flip(1:trainData.bin.dir.Nbin),  P.theta(1,1,b) = numel( trainData.bin.dir.ind{b} )/trainData.N;  end
for roi = 1:Nroi
    for z = 1:testData.N
        % Which activity bin did this unit/frame's activity fall into?
        %Ntemp = histcounts( log10(testData.activity(z,r)), Pdens.binEdges );
        %tempBinInd = find( Ntemp );
        tempBinInd = find( histcounts( log10(testData.activity(z,roi)), Pdens.binEdges ) );
        % How probable was that?
        P.a(z,roi) = Pdens.a{roi}(tempBinInd); %pdf( aPDF{r}, testData.activity(z,r) ); % diff( cdf( aPDF{r}, [testData.activity(z,r)-dA, testData.activity(z,r)+dA] ) );
        for b = 1:testData.bin.dir.Nbin
            P.a_theta(z,roi,b) = Pdens.a_theta{b,roi}(tempBinInd); % pdf( a_thetaPDF{b,r}, testData.activity(z,r) ); %diff( cdf( a_thetaPDF{b,r}, [testData.activity(z,r)-dA, testData.activity(z,r)+dA] ) ); % [0, testData.activity(z,r)+dA]
        end
    end
end
P.theta_a = (P.a_theta.*repmat(P.theta, testData.N, Nroi))./repmat(P.a, 1, 1, testData.bin.dir.Nbin); % P(thetaB|a) = P(a|thetaB)*P.theta(b)/P(a)  (assume units are independent)

% Check for pathological values
%fprintf('\n  P(theta) has %i pathological values', numel(find( isnan(P.theta) | isinf(P.theta))) );
%fprintf('\n  P(a) has %i pathological values', numel(find( isnan(P.a) | isinf(P.a))) ); %find( isnan(P.a) | isinf(P.a_theta) );
%fprintf('\n  P(a|theta) has %i pathological values', numel(find( isnan(P.a_theta) | isinf(P.a_theta))) ); %find( isnan(P.a_theta) | isinf(P.a_theta) );
%fprintf('\n  P(theta|a) has %i pathological values', numel(find( isnan(P.theta_a) | isinf(P.theta_a))) ); %find( isnan(P.theta_a) | isinf(P.theta_a) );

% Get the maximum likelihood theta
likeTheta = squeeze( sum( log10(P.theta_a), 2 ) ); % fprintf('\n  likeTheta has %i pathological values', numel(find( isnan(likeTheta) | isinf(likeTheta))) ); 
[~, maxLikeThetaBin] = max( likeTheta, [], 2 );
maxLikeThetaBinN = zeros(1, trainData.bin.dir.Nbin );
maxLikeThetaBinN(bUse) = histcounts( maxLikeThetaBin, Nuse);
maxLikeTheta = testData.bin.dir.cent(maxLikeThetaBin)';

% Assess overall accuracy/error
errStruct = struct('hit',[], 'Nhit',[], 'hitRate',NaN, 'HRbin',nan(1,trainData.bin.dir.Nbin), 'absErr',[], 'MAE',NaN, 'MAEbin',nan(1,trainData.bin.dir.Nbin), ...
    'nullHR',nan(Nboot,1), 'nullHRmean',NaN, 'nullHRci',nan(1,2), 'nullMAE',nan(Nboot,1), 'nullMAEmean',NaN, 'nullMAEci',nan(1,2), ...
    'nullHRbin',nan(Nboot,trainData.bin.dir.Nbin), 'nullMAEbin',nan(Nboot,trainData.bin.dir.Nbin)   );
errStruct.hit = maxLikeThetaBin - testData.trueBin == 0;
errStruct.Nhit = sum(errStruct.hit); %numel(errStruct.hit);
errStruct.hitRate = errStruct.Nhit/testData.N;
errStruct.absErr = abs(testData.dir - maxLikeTheta);
errStruct.absErr(errStruct.absErr > pi) = 2*pi - errStruct.absErr(errStruct.absErr > pi);
errStruct.absErr = rad2deg(abs(errStruct.absErr));
errStruct.MAE = median( errStruct.absErr );
fprintf('\n  Hit rate = %2.2f, median absolute error = %2.2f deg\n  ', errStruct.hitRate, errStruct.MAE )
% Assess binwise accuracy/error
for b = bUse
    binTestInd = find(testData.trueBin == b);
    errStruct.HRbin(b) = sum( maxLikeThetaBin(binTestInd) - testData.trueBin(binTestInd) == 0 )/numel(trueTestBinCell{b});
    errStruct.MAEbin(b) = median(errStruct.absErr(binTestInd)); %rad2deg(median( binAbsErr ));
end

% How accurate would random guessing be?
if Nboot > 0
    for n = 1:Nboot
        % Hit Rate
        nullBins = bUse( randi(Nuse, 1, testData.N) )'; % bUse
        nullHit = nullBins-testData.trueBin == 0;
        errStruct.nullHR(n) = sum(nullHit)/testData.N;
        % Absolute error
        nullAE = abs(testData.dir - testData.bin.dir.cent(nullBins)'); % nullBins*dTheta
        nullAE(nullAE > pi) = 2*pi - nullAE(nullAE > pi);
        nullAE = rad2deg(nullAE);
        errStruct.nullMAE(n) = median(nullAE);
        % Assess binwise accuracy/error
        for b = bUse
            binTestInd = find(testData.trueBin == b);
            errStruct.nullHRbin(n,b) = sum(nullHit(binTestInd))/numel(trueTestBinCell{b});
            errStruct.nullMAEbin(n,b) = median( nullAE(binTestInd) );
        end
    end
    errStruct.nullHRci = prctile(errStruct.nullHR, [5,95]);
    errStruct.nullHRmean = mean( errStruct.nullHR );
    errStruct.nullMAEci = prctile(errStruct.nullMAE, [5,95]);
    errStruct.nullMAEmean = mean( errStruct.nullMAE );
end
toc

if show
    FS = 14;
    MS = 4;
    opt1 = {[0.06,0.02], [0.09,0.04], [0.07, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
    opt2 = {[0.06,0.02], [0.01,0.04], [0.07, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
    close all; clearvars sp h;
    thetaHistLim = [testData.bin.dir.lim, pi]; %thetaHistLimTick = rad2deg(testData.bin.dir.cent(1:end-1))+180;
    thetaCentLimTick = unique( 180 + rad2deg( trainData.bin.dir.cent(1:end-1) ) );
    thetaBinLimTick = unique( 180 + rad2deg( trainData.bin.dir.lim ) );   %  
    MaxLikeBreakdown = figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'color','w' ); %  'WindowState','maximized'
    % Show the training/testing data
    %{
    subtightplot(3,3,[1,4,7],opt{:});
    if trainData.N + testData.N < 2000
        %polarplot( dirIn, magIn, 'k.', 'markerSize',MS ); hold on; % drawing too many points in polarplot crashses matlab
        polarplot( trainData.dir, trainData.mag, 'b.', 'markerSize',MS ); hold on; 
        polarplot( testData.dir, testData.mag, 'r.', 'markerSize',MS );
        thetalim( rad2deg(thetaLim) ) 
        legend( sprintf('Train (n = %i,  %i per bin)', trainData.N, trainData.N/Nuse), sprintf('Test (n = %i)', testData.N), 'Location','southOutside'); % sprintf('All (n = %i)', NinFrames), 
        %legend('All','Train','Test', 'Location','southOutside');
        rlim([0,1.5]);
        set(gca, 'thetatick', thetaBinLimTick, 'FontSize', FS );
        title( 'Deform Data, Binning, Restrictions and Subsets' );
    end
    %}
    % Show likelihood, predicted theta, and stretch mag of test data
    frameTicks = cellfun(@numel, testData.bin.dir.ind ); frameTicks = cumsum(frameTicks(frameTicks > 0));
    frameTicks = unique([1, frameTicks(1:end-1)]);
    
    sp(1) = subtightplot(3,2,1,opt1{:});
    plot( likeTheta(:,bUse), '.', 'markersize',MS); % (testDirSortInd)
    %lgd = legend(cellstr( num2str( testData.bin.dir.cent(1:end-1) ) )', 'Location',[0.65, 0.85, 0.055, 0.09]); % 
    ylabel('Likelihood of \theta');  title(sprintf('Hit Rate = %2.3f  (null = %2.3f)', errStruct.hitRate, mean(errStruct.nullHR)));
    set(gca, 'Xtick',frameTicks, 'XtickLabel',[], 'FontSize',FS, 'TickDir','out' )
    xtickangle(45)

    sp(2) = subtightplot(3,2,3,opt1{:});
    plot( [maxLikeTheta,testData.dir], '.', 'markersize',MS );  hold on;
    ylim(thetaLim);  %ylim([-3.2, 3.2]);
    lgd(1) = legend('Prediction', 'Data' , 'Location',[0.52, 0.61, 0.066, 0.045]); % , 'Location',[0.65, 0.85, 0.055, 0.09]  legend( 'Location','best');
    set(gca,'Ytick',testData.bin.dir.lim, 'FontSize',FS, 'Xtick',frameTicks, 'XtickLabel',[], 'TickDir','out'); % [-pi, -pi/2, 0, pi/2, pi]  'YtickLabel',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'}, 
    xtickangle(45)
    ylabel('\theta'); title(sprintf('Absolute Error (Med = %2.2f deg)', errStruct.MAE));
    
    sp(3) = subtightplot(3,2,5,opt1{:});
    plot( testData.mag, '.', 'markersize',MS ); hold on;
    ylabel('Deform Magnitude'); xlabel('Test Frame'); 
    title(sprintf('%i training samples (%i per bin), %i test samples (at least %i per bin)', trainData.N, trainData.N/Nuse, testData.N, min(cellfun(@numel, testData.bin.dir.ind)) ));
    set(gca, 'FontSize',FS, 'TickDir','out', 'Xtick',frameTicks); % , 'Xtick',[]
    xtickangle(45)
    linkaxes(sp,'x');
    xlim([-Inf,Inf]);

    % Distribution of predicted and observed directions
    subtightplot(3,2,2,opt2{:})
    polarplot( testData.bin.dir.cent(bUse), cellfun(@numel, testData.bin.dir.ind(bUse))/testData.N, 'k.' ); hold on;
    polarplot( testData.bin.dir.cent(bUse), maxLikeThetaBinN(bUse)/testData.N, 'ko' );
    %thetalim( rad2deg(thetaLim) );
    lgd(2) = legend('Test Data', 'Predicted', 'Location', [0.82, 0.9, 0.07, 0.04]); % [0.9, 0.9, 0.06, 0.04]
    title('Distribution of Stretches');
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', thetaCentLimTick, 'ThetaTickLabel',[] );
    
    %{
    polarhistogram( testData.dir,  thetaHistLim, 'normalization','probability' );
    thetalim( rad2deg(thetaLim) ) %thetalim( thetaLim )
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', thetaBinLimTick, 'RaxisLocation',-90  ); % thetaHistLimTick
    title('Distribution of Test Stretches');
    %}

    % Distribution of absolute error
    subtightplot(3,2,4,opt2{:})
    for b = bUse
        h(1) = polarplot( testData.bin.dir.cent(b)*[1,1], prctile( errStruct.nullHRbin(:,b), [5, 95] ), 'k' ); hold on;
        h(2) = polarplot( testData.bin.dir.cent(b), errStruct.HRbin(b), 'ko' );
    end
    lgd(3) = legend(h, {'5/95% Null CI', 'Prediction'}, 'Location', [0.82, 0.35, 0.08, 0.04]);
    %thetalim( rad2deg(thetaLim) );
    title('Binwise Hit Rates');
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', thetaCentLimTick, 'ThetaTickLabel',[] ); % , 'RaxisLocation',-90 
    
    subtightplot(3,2,6,opt2{:});%subplot(3,3,5);
    for b = bUse
        polarplot( testData.bin.dir.cent(b)*[1,1], prctile( errStruct.nullMAEbin(:,b), [5, 95] ), 'k' ); hold on;
        polarplot( testData.bin.dir.cent(b), errStruct.MAEbin(b), 'ko' );
    end
    %thetalim( rad2deg(thetaLim) );
    title('Binwise Median Absolute Error (deg)');
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', thetaCentLimTick, 'ThetaTickLabel',[] ); % , 'RaxisLocation',-90 

    if ~isempty(subName)
        if ~isinf(pauseTime),  pause(pauseTime);  else,  pause;  end
        figName = sprintf('MaxLikeBreakdown_%s.pdf', subName ); 
        figPath = ['C:\Users\ablaeser\Documents\Afferent Paper\MaxLikeResults\', figName];
        export_fig( figPath, '-pdf', '-painters','-q101', '-append', MaxLikeBreakdown );
        fprintf('\nSaved %s', figPath );
        %clf;
    end
end
end