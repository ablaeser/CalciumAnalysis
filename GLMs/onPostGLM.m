%% Use GLM to assess contribution of different variables
onPostResult = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\Loco\';
clearvars onPostPred onPostResp onPostResult onPostSumm
rLocoPostFit = cell(1,Nexpt);
for x = x3Dcsd %xPresent
    % GLMparallel options
    onPostOpts(x).show = true;
    onPostOpts(x).rShow = [];
    onPostOpts(x).figDir = figDir;
    onPostOpts(x).alpha = 0.01;  % The regularization parameter, default is 0.01
    onPostOpts(x).standardize = true; 
    onPostOpts(x).trainFrac = 0.75; % 1; %
    onPostOpts(x).Ncycle = 20;
    onPostOpts(x).distribution = 'gaussian'; % 'poisson';
    onPostOpts(x).CVfold = 10;
    onPostOpts(x).minDev = 0.1;  minDev = onPostOpts(x).minDev;
    onPostOpts(x).minDevFrac = 0.1;
    onPostOpts(x).maxP = 0.05;
    onPostOpts(x).Nshuff = 0;  Nshuff = onPostOpts(x).Nshuff;
    onPostOpts(x).minShuff = 15; 
    onPostOpts(x).window = [0, 0]; % [0,0]; % [-0.5, 0.5]; % 
    onPostOpts(x).lopo = true; %false; %
    onPostOpts(x).frameRate = expt(x).scanRate; 
    onPostOpts(x).minShuffFrame = round( onPostOpts(x).frameRate*onPostOpts(x).minShuff );
    windowFrame = round(onPostOpts(x).window*onPostOpts(x).frameRate);
    onPostOpts(x).shiftFrame = windowFrame(1):windowFrame(2);
    onPostOpts(x).maxShift = max( abs(windowFrame) );
    onPostOpts(x).Nshift = numel( onPostOpts(x).shiftFrame );  %Nshift = onPostOpts(x).Nshift;
    onPostOpts(x).lags = onPostOpts(x).shiftFrame/onPostOpts(x).frameRate;
    onPostOpts(x).name = strcat(expt(x).name, '_onPostGLM');

    % Concatenate pre-CSD bouts
    onSubCat = [];
    for v = 1:NallVars
        subData = onStruct(x).(allVars{v}).data(:,:,onStruct(x).csd.post.bout) - onStruct(x).(allVars{v}).pre(:,:,onStruct(x).csd.post.bout); % subtract each bout's pre-onset mean
        onSubCat.(allVars{v}) = reshape( permute( subData, [1,3,2]), ...
            size(onStruct(x).(allVars{v}).data,1)*onStruct(x).csd.post.Nbout, size(onStruct(x).(allVars{v}).data,2) );       
    end
    onSubCat.T = repmat( onStruct(x).T, onStruct(x).csd.post.Nbout, 1);
    
    % Define predictors
    onPostPred(x) = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % onPostPred(x).data = []; %onPostPred(x).data = onSubCat.T; onPostPred(x).name = {'Time'};
    for v = 2:NallVars
        onPostPred(x).data = [onPostPred(x).data, mean(onSubCat.(allVars{v}), 2, 'omitnan')]; %rand(1000,3); %
        onPostPred(x).name = [onPostPred(x).name, allVars(v)];
    end
    onPostPred(x).N = size(onPostPred(x).data,2);
    for p = flip(1:onPostPred(x).N), onPostPred(x).lopo.name{p} = ['No ',onPostPred(x).name{p}]; end
    onPostPred(x).fam.col = {1:2, 3:4, 5:6, 7:8, 9, 10};  %{1, 2:3, 4:5, 6:7, 8, 9}; 
    onPostPred(x).fam.N = numel(onPostPred(x).fam.col); 
    onPostPred(x).fam.name = {'Translation', 'Scale', 'Stretch', 'Shear', 'Z Shift', 'Speed'}; % 'Onset Time',
    
    % Define responses
    onPostResp(x).data = onSubCat.fluor; %[merged(x).scaleap-medAP, merged(x).scaleml-medML]; % , merged(x).dScAP, merged(x).dScML
    onPostResp(x).N = size(onPostResp(x).data,2); 
    onPostResp(x).name = repmat({'Fluor'}, 1, onPostResp(x).N ); %{'Scale-AP','Scale-ML'}; % , 'dScAP', 'dScML
    nanFrame = find(any(isnan([onPostPred(x).data, onPostResp(x).data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoved %i NaN-containing frames', numel(nanFrame));
    onPostPred(x).data(nanFrame,:) = []; onPostResp(x).data(nanFrame,:) = [];
    
    % Run the GLM
    onPostOpts(x).load = false; %  true; %
    onPostOpts(x).saveRoot = sprintf('%s%s', expt(x).dir, onPostOpts(x).name  ); %''; %
    [onPostResult{x}, onPostSumm(x), ~, onPostPred(x), onPostResp(x)] = GLMparallel(onPostPred(x), onPostResp(x), onPostOpts(x)); % onPostResult{x}, onPostSumm(x), ~, onPostPred(x), onPostResp(x)
    rLocoPostFit{x} = intersect(onStruct(x).fluor.responder, onPostSumm(x).rFit);
    nLocoPreFit(x) = numel( rLocoPostFit{x} );
end

%% Summarize GLM results
close all; clearvars sp SP;
OnPostGLMresults = figure('WindowState','maximized', 'color','w');
leftOpt = {[0.02,0.07], [0.06,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
rightOpt = {[0.1,0.07], [0.1,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
jitterWidth = 0.45;
xAngle = 30;
Ncol = 3;
for x = xCSD %xPresent
    Nrow = onPostPred(x).N+1; 
    spGrid = reshape( 1:Nrow*Ncol, Ncol, Nrow )';
    subtightplot(onPostPred(x).N+1, 3, 1:2, leftOpt{:});
    imagesc( onPostResp(x).data' );
    ylabel('Fluor', 'Interpreter','none');
    set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'Ytick',onStruct(x).fluor.responder, 'Xtick',[]);
    text( repmat(size(onPostResp(x).data,1)+1, onPostSumm(x).Nfit, 1), onPostSumm(x).rFit+0.5, '*', 'VerticalAlignment','middle', 'FontSize',8);
    impixelinfo;
    
    for v = 1:onPostPred(x).N
        subtightplot(Nrow, Ncol, spGrid(v+1,1:2), leftOpt{:});
        plot( onPostPred(x).data(:,v) ); hold on;
        ylabel(onPostPred(x).name{v}, 'Interpreter','none');
        xlim([-Inf,Inf]);
        if v < onPostPred(x).N
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'XtickLabel',[]);
        else
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off');
        end
    end
    xlabel('Scan');
    
    subtightplot(3,3,3, rightOpt{:});
    bar([numel(onStruct(x).fluor.responder), onPostSumm(x).Nfit, numel(rLocoPostFit{x})]/expt(x).Nroi );
    set(gca,'Xtick',1:3, 'XtickLabel',{'Loco','Fit','Both'}, 'box','off');
    ylabel('Fraction of ROI');
    ylim([0,1]);
    
    subtightplot(3,3,6, rightOpt{:});
    JitterPlot( onPostSumm(x).lopo.devFrac(onPostSumm(x).rFit,:), jitterWidth ); hold on;
    line([0,onPostPred(x).N], [1,1], 'color','k', 'lineStyle','--');
    xlim([0,onPostPred(x).N]);
    ylabel('Fraction of total deviance'); title('Leave One Predictor Out (well-fit units only)');
    set(gca, 'Xtick',1:onPostPred(x).N,  'XtickLabel', onPostSumm(x).lopo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    subtightplot(3,3,9, rightOpt{:});
    JitterPlot( onPostSumm(x).lofo.devFrac(onPostSumm(x).rFit,:), jitterWidth ); hold on;
    line([0,onPostPred(x).fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
    xlim([0,onPostPred(x).fam.N]+0.5);
    ylabel('Fraction of total deviance'); title('Leave One Family Out (well-fit units only)');
    set(gca, 'Xtick',1:onPostPred(x).fam.N,  'XtickLabel', onPostSumm(x).lofo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    % {
    figPath = sprintf('%s%s.tif', figDir, onPostOpts(x).name);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    print(OnPostGLMresults, figPath, '-dtiff' ); pause(1);   
    clf;
    %}
    %pause; clf;
end


%% Plot coefficient values for each predictor, unit and experiment
close all; clearvars sp SP;
OnPostGLMcoeff = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.04], [0.09,0.05], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
for v = 1:onPostPred(x).N
    c = 0;
    for x = find(~cellfun(@isempty, rLocoPostFit)) %xPresent
        c = c+1;
        subtightplot(1, onPostPred(x).N, v, opt{:});
        plot(c, onPostSumm(x).peakCoeff(rLocoPostFit{x},v), 'k.' ); hold on;
    end
    line([0,c+1], [0,0], 'color','k');
    title( onPostPred(x).name{v}, 'Interpreter','none' );
    if v == 1, ylabel('Coefficient'); end
    xlim([0,c+1]);
    set(gca, 'Xtick', 1:c, 'XtickLabel', {expt(xPresent).name}, 'TickLabelInterpreter','none' );
    xtickangle(45);
    axis square;
end

%% LOFO/LOPO results, pooling well-fit, locomotion-sensitive ROI across experiments
lopoPool = []; lofoPool = [];
for x = find(~cellfun(@isempty, rLocoPostFit))
    lopoPool = vertcat( lopoPool, onPostSumm(x).lopo.devFrac(rLocoPostFit{x}, :) );
    lofoPool = vertcat( lofoPool, onPostSumm(x).lofo.devFrac(rLocoPostFit{x}, :) );
end
close all; clearvars sp SP;
OnPostLeaveOneOutResults = figure('WindowState','maximized', 'color','w');
opt = {[0.03,0.07], [0.09,0.05], [0.08,0.07]};  % {[vert, horz], [bottom, top], [left, right] }

subtightplot(1,2,1,opt{:});
JitterPlot( lopoPool, 0.45 );
line([0,onPostPred(xCSD(1)).N]+0.5, [1,1], 'color','k', 'lineStyle','--');
xlim([0,onPostPred(xCSD(1)).N]+0.5);
xtickangle(xAngle);
set(gca,'Xtick',1:onPostPred(xCSD(1)).N, 'XtickLabel',onPostPred(xCSD(1)).name, 'TickLabelInterpreter','none' );
title('Leave One Predictor Out (locomotion-sensitive, well-fit units only)'); ylabel('Fraction of total deviance');
axis square;

subtightplot(1,2,2,opt{:});
JitterPlot( lofoPool, 0.45 );
line([0,onPostPred(xCSD(1)).fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
xlim([0,onPostPred(xCSD(1)).fam.N]+0.5);
set(gca,'Xtick',1:onPostPred(xCSD(1)).fam.N, 'XtickLabel',onPostPred(xCSD(1)).fam.name, 'TickLabelInterpreter','none' );
title('Leave One Family Out'); ylabel('Fraction of total deviance');
axis square;

figPath = sprintf('%sLeaveOneOutResults_onsetPreCSD.tif', figDir);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(OnPostLeaveOneOutResults, figPath, '-dtiff' ); %pause(1); clf;