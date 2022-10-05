function [ stimGroup, respGroup, fitResult, fitSumm, sigmFit ] = StimResponse( stim, response, varargin )
%StimResponse 
Nresponse = size(response,2);
IP = inputParser;
addRequired( IP, 'stim', @isnumeric )
addRequired( IP, 'response', @isnumeric )
addOptional( IP, 'frameDelay', 0, @isnumeric )
addOptional( IP, 'groupPct', 10, @isnumeric )
addOptional( IP, 'Nsigma', 0, @isnumeric )
addParameter( IP, 'alpha', 0.05, @isnumeric )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'fit', true, @islogical )
addParameter( IP, 'setr', 1:Nresponse, @isnumeric )
addParameter( IP, 'names', {'Stimulus','Response'}, @iscell )
parse( IP, stim, response, varargin{:} ); 
frameDelay = IP.Results.frameDelay;
groupPct = IP.Results.groupPct;
Nsigma = IP.Results.Nsigma;
alphaVal = IP.Results.alpha;
show = IP.Results.show;
setr = IP.Results.setr;
fitSigmoid = IP.Results.fit;
varNames = IP.Results.names;
% setr = 1:Nresponse; Nsigma = 0;

if groupPct > 0 && Nsigma > 0, error('Nsigma and groupPct cannot both be non-zero');  end

% Suppress missing data
missingStimFrame = find(isnan(sum(stim,2)));
missingRespFrame = find(isnan(sum(response,2)));
suppressFrame =  unique(vertcat(missingStimFrame, missingStimFrame + frameDelay, missingRespFrame, missingRespFrame-frameDelay) );
stim(suppressFrame,:) = [];
response(suppressFrame,:) = [];

% Group data by stimulus intensity
if Nsigma > 0
    stimGroup.type = 'Sigma';
    stimMean = nanmean(stim); stimStd = nanstd(stim);
    Ngroup = 2*Nsigma+1;
    for g = -Nsigma:Nsigma
        stimGroup.lim(g+Nsigma+1,:) = [(stimMean+(g-1)*(stimStd)), (stimMean+g*stimStd)];
        stimGroup.ind{g+Nsigma+1,1} = find( stim >= (stimMean+(g-1)*(stimStd)) & stim < (stimMean+g*stimStd) ); 
    end
elseif groupPct > 0
    stimGroup.type = 'Percentile';
    Ngroup = ceil(100/groupPct); 
    for g = flip(1:Ngroup)
        stimGroup.lim(g,:) = [prctile(stim,(g-1)*groupPct), prctile(stim,g*groupPct)];
        stimGroup.ind{g,1} = find( stim >= prctile(stim,groupPct*(g-1)) & stim < prctile(stim,groupPct*g) ); 
    end
else
    error('Grouping not defined');
end  

% Group responses with delay
for g = flip(1:Ngroup)
    respGroup.ind{g,1} = stimGroup.ind{g} + frameDelay; 
    respGroup.ind{g,1}(respGroup.ind{g} > size(stim,1) | respGroup.ind{g} < 1) = []; %exclude indices that are shifted out of the dataset
    respGroup.val{g,1} = response(respGroup.ind{g},:);
    %respGroup.mean(g,:) = mean( respGroup.val{g}, 1 );
    %respGroup.std(g,:) = std( respGroup.val{g}, 0, 1 );
    %respGroup.sem(g,:) = respGroup.std(g,:)/sqrt(respGroup.N(g));
    stimGroup.ind{g,:} = respGroup.ind{g} - frameDelay; %exclude indices that are shifted out of the dataset
    stimGroup.val{g,:} = stim(stimGroup.ind{g});
    %stimGroup.mean(g,:) = mean( stimGroup.val{g} );
end
respGroup.N = cellfun(@numel, respGroup.ind );

% Remove empty groups
emptyGroup = find(respGroup.N == 0); %find(stimGroup.N == 0);  % 
if ~isempty(emptyGroup)
    fprintf('\nWarning: Grouping includes %i empty groups!\n', sum(respGroup.N == 0) ); 
    stimGroup.lim(emptyGroup,:) = [];
    stimGroup.ind(emptyGroup,:) = [];
    stimGroup.val(emptyGroup,:) = [];
    %stimGroup.mean(emptyGroup,:) = [];

    respGroup.ind(emptyGroup,:) = [];
    respGroup.N(emptyGroup) = [];
    respGroup.val(emptyGroup,:) = [];
    %respGroup.mean(emptyGroup,:) = [];
    %respGroup.std(emptyGroup,:) = [];
    %respGroup.sem(emptyGroup,:) = [];
end
Ngroup = numel(respGroup.N);
stimGroup.mean = cellfun(@mean, stimGroup.val);
respGroup.mean = cellfun(@mean, respGroup.val, 'UniformOutput',false);
respGroup.mean = vertcat(respGroup.mean{:});
respGroup.std = cellfun(@std, respGroup.val, 'UniformOutput',false);
respGroup.std = vertcat(respGroup.std{:});
respGroup.sem = respGroup.std./sqrt(respGroup.N);

% what's the lowest group that's significantly different from baseline (Multiple comparisons test)
respGroup.gFirstSig = nan(1,Nresponse); respGroup.minStim = nan(1,Nresponse);
for r = setr 
    % Gather all the responses from all groups, for a given unit
    tempVals = nan( max(respGroup.N), Ngroup );
    for g = flip(1:Ngroup), tempVals(1:respGroup.N(g),g) = respGroup.val{g}(:,r); end
    % Are any groups significantly different?
    [pGroupDiff, ~, GroupDiffStats] = anova1( tempVals, [], 'off'  ); %kruskalwallis( cell2padmat(respGroup.val), [], 'off' ); %  GroupDiffTbl
    if pGroupDiff < alphaVal % Find the minimal group that's significantly different from the first group
        MCtbl = multcompare( GroupDiffStats, 'CType','dunn-sidak', 'display','off' ); %
        MCtbl = MCtbl(1:Ngroup-1,:);
        tempMinStimInd = find( MCtbl(:,6) < alphaVal, 1, 'first' );
        if ~isempty(tempMinStimInd)
            respGroup.gFirstSig(r) = MCtbl(tempMinStimInd,2);
            respGroup.minStim(r) = stimGroup.mean(respGroup.gFirstSig(r));
        end
    end
end

% Fit to sigmoid (optional)
sigmFit = cell(Nresponse,1); 
fitResult = repmat( struct('R2',NaN, 'thresh',NaN, 'pThresh',NaN, 'rate',NaN, 'pRate',NaN, 'factor',NaN, 'pFactor',NaN, 'offset',NaN, 'pOffset',NaN ), 1, Nresponse );
fitSumm = struct('rGood',[], 'Ngood',NaN, 'goodFrac',NaN, 'R2good',[], 'threshGood',[], 'rateGood',[]);
if fitSigmoid || nargout > 4 % true %
    modelSigm = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4) ); % beta = [factor, growth rate, threshold, constant response offset]
    %tic
    for r = setr 
        % Determine initial parameters for sigmoidal fit
        % Constant offset response
        betaInit(4) = respGroup.mean(1,r); %min( respGroup.mean(:,r) ); 
        % Saturation response
        betaInit(1) = respGroup.mean(end,r) - betaInit(4); % max(respGroup.mean(:,r))- betaInit(4); 
        % Response growth rate
        [~,minInd] = min( abs(respGroup.mean(:,r) - (betaInit(1)/2 + betaInit(4))) );
        betaInit(2) = stimGroup.mean(minInd); 
        % Threshold stimulus
        if minInd > 1 && minInd < numel(stimGroup.mean)
            betaInit(3) = stimGroup.mean(minInd+1) - stimGroup.mean(minInd-1);
        elseif minInd == 1
            betaInit(3) = stimGroup.mean(minInd+1) - stimGroup.mean(minInd);
        elseif minInd == numel(stimGroup.mean)
            betaInit(3) = stimGroup.mean(minInd) - stimGroup.mean(minInd-1); 
        else 
            betaInit(3) = 1;
        end
        sigmFit{r} = fitnlm( stimGroup.mean, respGroup.mean(:,r), modelSigm, betaInit );
        fitResult(r).R2 =  sigmFit{r}.Rsquared.Adjusted; % goodness of fit
        fitResult(r).thresh = sigmFit{r}.Coefficients.Estimate(3);
        fitResult(r).pThresh = sigmFit{r}.Coefficients.pValue(3);
        fitResult(r).rate = sigmFit{r}.Coefficients.Estimate(2);
        fitResult(r).pRate = sigmFit{r}.Coefficients.pValue(2);
        fitResult(r).factor = sigmFit{r}.Coefficients.Estimate(1);
        fitResult(r).pFactor = sigmFit{r}.Coefficients.pValue(1);
        fitResult(r).offset = sigmFit{r}.Coefficients.Estimate(4);
        fitResult(r).pOffset = sigmFit{r}.Coefficients.pValue(4);
    end
    %toc
    goodR2 = 0.9; goodP = 0.05;
    R2pool = [fitResult.R2]; threshPool = [fitResult.thresh]; ratePool = [fitResult.rate];
    % Summarize results of sigmoidal fitting
    fitSumm.rGood = find( R2pool > goodR2 & threshPool > 0 & threshPool < stimGroup.mean(end) & ratePool > 0 & [fitResult.factor] > 0 & [fitResult.offset] > 0 & [fitResult.pThresh] < goodP & [fitResult.pRate] < goodP );
    fitSumm.Ngood = numel(fitSumm.rGood);
    fitSumm.goodFrac = fitSumm.Ngood/numel(setr);
    fitSumm.R2good = R2pool( fitSumm.rGood );
    fitSumm.threshGood = threshPool( fitSumm.rGood );
    fitSumm.rateGood = ratePool( fitSumm.rGood );
    %fprintf('\nGoodness-of-fit criteria: R2adj > %1.2f;  pThresh & pRate < %1.2f;  all params > 0;  threshold < max stimulus \nFound %i of %i responses (%2.1f pct) meeting these criteria\n', ...
    %    goodR2, goodP, fitSumm.Ngood, numel(setr), 100*fitSumm.goodFrac );
end

% Plot the process (optional)
if show
    figure('WindowState', 'maximized', 'color','w', 'PaperOrientation','landscape');
    subplot(2,2,1);
    plot( stim ); hold all;
    for g = 1:Ngroup
        plot( stimGroup.ind{g}, stim(stimGroup.ind{g}), '.' );
        line([1, numel(stim)], stimGroup.lim(g,1)*[1,1], 'LineStyle','--', 'Color','k');
        line([1, numel(stim)], stimGroup.lim(g,2)*[1,1], 'Color','k');
    end
    axis tight;
    xlabel('Frame'); ylabel(varNames{1}); title(sprintf('%s (%i)', stimGroup.type, groupPct+Nsigma));
    
    for r = fitSumm.rGood %setr %goodResp % 1:Nresponse
        subplot(2,2,3); cla;
        plot( response(:,r) ); hold all;
        for g = 1:Ngroup,  plot( respGroup.ind{g}, respGroup.val{g}(:,r) , '.' );  end % response(respGroup.ind{g},r)
        xlabel('Frame'); ylabel( varNames{2} ); title(sprintf('r = %i', r )); axis tight; 
        
        subplot(2,2,[2,4]); cla;
        errorbar( stimGroup.mean, respGroup.mean(:,r), respGroup.sem(:,r), 'LineStyle','none', 'Color','k' ); hold on;
        if ~isnan( respGroup.gFirstSig(r) )
            errorbar( stimGroup.mean(respGroup.gFirstSig(r)), respGroup.mean(respGroup.gFirstSig(r),r), respGroup.sem(respGroup.gFirstSig(r),r), 'color','b', 'LineStyle','none' );
        end
        if fitSigmoid || nargout > 4
            xRange = linspace(0,stimGroup.mean(end))';  %(0:0.01:stimGroup.mean(end))';
            plot( xRange, predict(sigmFit{r}, xRange), 'k' ); 
            title( sprintf('R^2adj = %2.3f, Threshold = %2.5f (p=%1.3f), Rate = %2.5f (p=%1.3f)', fitResult(r).R2, fitResult(r).thresh, fitResult(r).pThresh, fitResult(r).rate, fitResult(r).pRate ) );
        end
        xlabel( sprintf('%s (%2.1f %s grouping)', varNames{1}, groupPct, stimGroup.type )); ylabel( sprintf('%s (%i frame delay)', varNames{2}, frameDelay ) );
        axis square;
        if r < Nresponse, pause; end
    end
end
end