%% Use GLM to assess contribution of different variables
preCSD_pred = cell(1,Nexpt); preCSD_resp = cell(1,Nexpt); preCSD_opts = cell(1,Nexpt); preCSD_result = cell(1,Nexpt); preCSD_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
%rLocoPreFit = cell(1,Nexpt);
GLMname = 'postCSD_GLM_3D';
for x = x3Dcsd % xPresent %
    % GLMparallel options
    postCSD_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_postCSDglm');
    %postCSD_opts{x}.show = true;
    postCSD_opts{x}.rShow = NaN;
    postCSD_opts{x}.figDir = ''; % figDir;
    postCSD_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    postCSD_opts{x}.standardize = true; 
    postCSD_opts{x}.trainFrac = 0.75; % 1; %
    postCSD_opts{x}.Ncycle = 20;
    postCSD_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    postCSD_opts{x}.CVfold = 10;
    postCSD.opts(x).nlamda = 1000;
    postCSD.opts(x).maxit = 5*10^5;
    postCSD_opts{x}.minDev = 0.1;  minDev = postCSD_opts{x}.minDev;
    postCSD_opts{x}.minDevFrac = 0.1;
    postCSD_opts{x}.maxP = 0.05;
    postCSD_opts{x}.Nshuff = 0;  Nshuff = postCSD_opts{x}.Nshuff;
    postCSD_opts{x}.minShuff = 15; 
    postCSD_opts{x}.window = [-5,5]; % [0,0]; % [-0.5, 0.5]; % 
    postCSD_opts{x}.lopo = true; %false; %
    postCSD_opts{x}.frameRate = expt(x).scanRate; 
    postCSD_opts{x}.minShuffFrame = round( postCSD_opts{x}.frameRate*postCSD_opts{x}.minShuff );
    windowFrame = round(postCSD_opts{x}.window*postCSD_opts{x}.frameRate);
    postCSD_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    postCSD_opts{x}.maxShift = max( abs(windowFrame) );
    postCSD_opts{x}.Nshift = numel( postCSD_opts{x}.shiftFrame );  %Nshift = postCSDOpts(x).Nshift;
    postCSD_opts{x}.lags = postCSD_opts{x}.shiftFrame/postCSD_opts{x}.frameRate;

    % Concatenate input variables pre-CSD
    if ~isnan(expt(x).csd)
        postCSDruns = expt(x).csd:expt(x).Nruns;
        tempStateCat = vertcat(loco{x}(postCSDruns).stateDown);
        tempTransVel = abs( mean( vertcat(deform{x}(postCSDruns).DtransMag), 2, 'omitnan') );
        tempScaleMag = mean( vertcat(deform{x}(postCSDruns).scaleMag), 2, 'omitnan');
        tempStretchMag = mean( vertcat(deform{x}(postCSDruns).stretchMag), 2, 'omitnan');
        tempShearMag = mean( vertcat(deform{x}(postCSDruns).shearMag), 2, 'omitnan');
        tempDshearMag = mean( vertcat(deform{x}(postCSDruns).DshearMag), 2, 'omitnan');
        tempShift = vertcat(deform{x}(postCSDruns).shiftZ);
        tempShiftMean =  mean( tempShift(:,4:end-3), 2, 'omitnan' );
        tempCompress = vertcat(deform{x}(postCSDruns).compressZ);
        tempCompressMean = mean( tempCompress(:,4:end-2), 2, 'omitnan' );
        tempFluor = [fluor{x}.z];
        tempFluorCat = vertcat(tempFluor(postCSDruns).ROI);

        % Define predictors
        postCSD_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
        postCSD_pred{x}.data = [tempTransVel, tempStretchMag, tempShiftMean, tempStateCat]; %
        postCSD_pred{x}.name = {'TransSpeed', '|Stretch|', 'Z Shift', 'Locomotive State'}; % ,'Speed',
        postCSD_pred{x}.N = size(postCSD_pred{x}.data,2);
        for p = flip(1:postCSD_pred{x}(r).N), postCSD_pred{x}(r).lopo.name{p} = ['No ',postCSD_pred{x}(r).name{p}]; end
        % Set up leave-one-family-out
        postCSD_pred{x}.fam.col = {1:3, 4}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
        postCSD_pred{x}.fam.N = numel(postCSD_pred{x}.fam.col); 
        postCSD_pred{x}.fam.name = {'No Deformation', 'No Loco State'}; %{'All'};%  'Onset Time',

        % Define response
        postCSD_resp{x}.data = tempFluorCat; % tempScaleMag; %
        postCSD_resp{x}.N = size(postCSD_resp{x}.data,2); 
        postCSD_resp{x}.name = sprintfc('Fluor %i', 1:postCSD_resp{x}.N);

        % Remove scans with missing data 
        nanFrame = find(any(isnan([postCSD_pred{x}.data, postCSD_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
        fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
        postCSD_pred{x}.data(nanFrame,:) = []; postCSD_resp{x}.data(nanFrame,:) = [];

        % Run the GLM
        postCSD_opts{x}.load = true; % false; %  
        postCSD_opts{x}.saveRoot = sprintf('%s%s', expt(x).dir, postCSD_opts{x}.name  ); %''; %

        [postCSD_result{x}, postCSD_summary{x}, ~, postCSD_pred{x}, postCSD_resp{x}] = GLMparallel(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}); % postCSD_result{x}, postCSD_summary(x), ~, postCSDPred(x), postCSDResp(x)
        %GLMresultFig = ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
        %pause;
        %rLocoPreFit{x} = intersect(onStruct(x).fluor.responder, postCSD_summary{x}.rGood);
        %nLocoPreFit(x) = numel( rLocoPreFit{x} );
    end
    
end

%% Compare GLM to data for each experiment
close all; clearvars sp SP;
PreGLMresults = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.06,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
rightOpt = {[0.1,0.07], [0.1,0.03], [0.04,0.02]};  % {[vert, horz], [bottom, top], [left, right] }\
jitterWidth = 0.45;
xAngle = 30;
Nrow = postCSD_pred{x}(r).N+1; Ncol = 3;
spGrid = reshape( 1:Nrow*Ncol, Ncol, Nrow )';
for x = find(~cellfun(@isempty,postCSD_pred))%xPresent
    subtightplot(postCSD_pred{x}(r).N+1, 3, 1:2, opt{:});
    imagesc( postCSD_resp{x}.data' );
    ylabel('Fluor', 'Interpreter','none');
    title( sprintf('%s', expt(x).name), 'Interpreter','none');
    set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'Xtick',[]); % , 'Ytick',onStruct(x).fluor.responder
    text( repmat(size(postCSD_resp{x}(r).data,1)+1, postCSD_summary{x}.Ngood, 1), postCSD_summary{x}.rGood+0.5, '*', 'VerticalAlignment','middle', 'FontSize',8);
    impixelinfo;
    
    for v = 1:postCSD_pred{x}(r).N
        subtightplot(Nrow, Ncol, spGrid(v+1,1:2), opt{:});
        plot( postCSD_pred{x}(r).data(:,v) ); hold on;
        ylabel(postCSD_pred{x}(r).name{v}, 'Interpreter','none');
        xlim([-Inf,Inf]);
        if v < postCSD_pred{x}(r).N
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off', 'XtickLabel',[]);
        else
            set(gca,'TickDir','out', 'TickLength',[0.003,0], 'box','off');
        end
    end
    xlabel('Scan');
    
    subtightplot(3,3,3, rightOpt{:});
    bar([postCSD_summary{x}.Ngood]/expt(x).Nroi ); % numel(onStruct(x).fluor.responder),   , numel(rLocoPreFit{x})
    set(gca,'Xtick',1, 'XtickLabel',{'Fit'}, 'box','off'); % 'Loco','Fit','Both'  :3
    ylabel('Fraction of ROI');
    ylim([0,1]);
    
    subtightplot(3,3,6, rightOpt{:});
    JitterPlot( postCSD_summary{x}.lopo.devFrac(:,postCSD_summary{x}.rGood)', jitterWidth ); hold on;
    line([0,postCSD_pred{x}(r).N+1], [1,1], 'color','k', 'lineStyle','--');
    xlim([0,postCSD_pred{x}(r).N+1]); ylim([0,Inf]); 
    ylabel('Fraction of total deviance'); title('Leave One Predictor Out (well-fit units only)');
    set(gca, 'Xtick',1:postCSD_pred{x}(r).N,  'XtickLabel', postCSD_summary{x}.lopo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    
    
    subtightplot(3,3,9, rightOpt{:});
    JitterPlot( postCSD_summary{x}.lofo.devFrac(:,postCSD_summary{x}.rGood)', jitterWidth ); hold on;
    line([0,postCSD_pred{x}(r).fam.N]+0.5, [1,1], 'color','k', 'lineStyle','--');
    xlim([0,postCSD_pred{x}(r).fam.N]+0.5);
    ylabel('Fraction of total deviance'); title('Leave One Family Out (well-fit units only)');
    set(gca, 'Xtick',1:postCSD_pred{x}(r).fam.N,  'XtickLabel', postCSD_summary{x}.lofo.name, 'TickDir','out', 'TickLength',[0.003,0], 'TickLabelInterpreter','none', 'box','off' ); 
    xtickangle(xAngle);
    ylim([0,Inf]);
    
    % {
    figPath = sprintf('%s%s_Deviance.tif', figDir, GLMname);
    if exist(figPath,'file'), delete(figPath); end
    fprintf('\nSaving %s', figPath);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', LocoSensitivePrePost ); pause(1);
    %print(PreGLMresults, figPath, '-dtiff' ); 
    pause%(1);   
    clf;
    %}
    %pause; clf;
end

%% Divide units into insensitive, mixed, loco-only or deformation-only units
postCSD_Nsubtype = []; k = 0;
for x = find(~cellfun(@isempty,postCSD_summary)) %xPresent
    postCSD_summary{x}.rIns = setdiff(1:expt(x).Nroi, postCSD_summary{x}.rGood);
    postCSD_summary{x}.nIns = numel( postCSD_summary{x}.rIns );
       
    postCSD_summary{x}.rDeform = setdiff(postCSD_summary{x}.lofo.rDependent{1}, postCSD_summary{x}.lofo.rDependent{2});
    postCSD_summary{x}.nDeform = numel( postCSD_summary{x}.rDeform );
    
    postCSD_summary{x}.rLoco = setdiff(postCSD_summary{x}.lofo.rDependent{2}, postCSD_summary{x}.lofo.rDependent{1});
    postCSD_summary{x}.nLoco = numel( postCSD_summary{x}.rLoco );
    
    postCSD_summary{x}.rMixed = setdiff(1:expt(x).Nroi, [postCSD_summary{x}.rIns, postCSD_summary{x}.rDeform, postCSD_summary{x}.rLoco]  );
    %postCSD_summary{x}.rMixed = intersect(postCSD_summary{x}.lofo.rDependent{1}, postCSD_summary{x}.lofo.rDependent{2});
    postCSD_summary{x}.nMixed = numel( postCSD_summary{x}.rMixed );
    
    k = k+1;
    postCSD_Nsubtype(k,:) = [postCSD_summary{x}.nIns, postCSD_summary{x}.nMixed, postCSD_summary{x}.nDeform, postCSD_summary{x}.nLoco]; % /expt(x).Nroi
end
postCSD_Nsubtype(k+1,:) = sum(postCSD_Nsubtype, 1);
postCSD_subtypeFrac = postCSD_Nsubtype./repmat( sum(postCSD_Nsubtype,2), 1, 4);

bar(postCSD_subtypeFrac,'stacked')

%% Show single examples of each subtype
for x = find(~cellfun(@isempty,postCSD_summary)) % xPresent
    postCSD_opts{x}.rShow = postCSD_summary{x}.rMixed;
    ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
    
    %postCSD_opts{x}.rShow = postCSD_summary{x}.rDeform;
    %ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
    
    %postCSD_opts{x}.rShow = postCSD_summary{x}.rLoco;
    %ViewGLM(postCSD_pred{x}, postCSD_resp{x}, postCSD_opts{x}, postCSD_result{x}, postCSD_summary{x});
end

%% Pool results across experiments
postCSDdevPool = []; preCSDdevPool = []; 
preCSDcoeffPool = []; postCSDcoeffPool = []; 
postGoodDevPool = [];
postCSDdevFracPool = []; %lofoDevFracPool = [];
for x = find(~cellfun(@isempty,postCSD_summary))%xPresent%find(~cellfun(@isempty, rLocoPreFit))
    %if ~isempty( postCSD_summary{x}.rGood )
        preCSDdevPool = [preCSDdevPool, preCSD_summary{x}.dev]; 
        postCSDdevPool = [postCSDdevPool, postCSD_summary{x}.dev]; 
        postGoodDevPool = [postGoodDevPool, postCSD_summary{x}.dev( postCSD_summary{x}.rGood )];
        postCSDdevFracPool = [postCSDdevFracPool, vertcat(postCSD_summary{x}.lopo.devFrac(:,postCSD_summary{x}.rGood), postCSD_summary{x}.lofo.devFrac(:,postCSD_summary{x}.rGood) )]; % rLocoPreFit{x}, :
    
        preCSDcoeffPool = vertcat(preCSDcoeffPool, preCSD_summary{x}.peakCoeff );
        postCSDcoeffPool = vertcat(postCSDcoeffPool, postCSD_summary{x}.peakCoeff );
        %end
end
%% Summarize deviance explained
opt = {[0.02,0.07], [0.1,0.07], [0.09,0.09]};  % {[vert, horz], [bottom, top], [left, right] }
DevianceFig = figure('WindowState','maximized', 'color','w');
k = 1; clearvars h;
exptColor = distinguishable_colors(numel(xPresent));
subtightplot(1,3,1,opt{:});
for x = find(~cellfun(@isempty,postCSD_summary))% xPresent
    [Ftemp, Xtemp] = ecdf( postCSD_summary{x}.dev ); hold on;
    h(k) = plot(Xtemp, Ftemp, 'color',exptColor(k,:) );
    k = k + 1;
end
[Fdev, Xdev] = ecdf( postCSDdevPool );
h(k) = plot( Xdev, Fdev, 'color','k', 'LineWidth',2 ); 
axis square;
legend(h, {expt(xPresent).name, 'Pooled', 'Threshold'}, 'Location','southEast', 'Interpreter','none', 'AutoUpdate',false );
xlim([0, 0.6]);
line(postCSD_opts{xPresent(1)}.minDev*[1,1], [0,1], 'Color','r', 'LineStyle','--'); % h(k+1) = 
xlabel('Deviance Explained'); ylabel('Fraction of Units');
title( sprintf('%s Fit Results', GLMname), 'Interpreter','none' );

subtightplot(1,3,2,opt{:});
JitterPlot( 1 - postCSDdevFracPool', 0.45, 'ErrorCap',10, 'monochrome',0.6); hold on;
line([0,size(postCSDdevFracPool,1)+1], [0,0], 'color','k');
axis square;
ylabel('Relative Explanatory Value'); %ylabel('Cumulative Distribution');
ylim([-1,1]);
set(gca,'Xtick', 1:size(postCSDdevFracPool,1), 'XtickLabel',[postCSD_pred{x}.name, postCSD_pred{x}.fam.name]);
xtickangle(30);
title('Well-Fit Units');

subtightplot(1,3,3,opt{:});
bar(postCSD_subtypeFrac,'stacked');
ylabel('Fraction of Units');
barPos = get(gca, 'Position');
xlim([0, size(postCSD_subtypeFrac,1)+1]);
axis square;
title('Subtype Breakdown');
legend('Insensitive','Mixed','Deform-only','Loco-only', 'Location','EastOutside');
set(gca, 'Xtick',1:size(postCSD_subtypeFrac,1), 'XtickLabel',{expt(xPresent).name, 'Pooled'}, 'TickLabelInterpreter','none', 'FontSize',10, 'Position',barPos );
xtickangle(30);

% Save the figure
figPath = sprintf('%s%sDevianceResults.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
fprintf('\nSaving %s', figPath);
print(DevianceFig, figPath, '-dtiff', '-r300'); %pause(1); clf;

%% Show mean coeff by type
Ntype = 4;
typeColor = distinguishable_colors(Ntype);
close all;
SubtypeCoeffFig = figure('WindowState','maximized', 'color','w');
opt = {[0.02,0.07], [0.1,0.1], [0.1,0.06]};  % {[vert, horz], [bottom, top], [left, right] }
LW = 1.5;
colororder( typeColor ) 
for x = find(~cellfun(@isempty,postCSD_summary))%xPresent
    tempInsCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rIns ).coeff );
    meanInsCoeff = mean(tempInsCoeff, 3, 'omitnan' );
    
    tempMixedCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rMixed ).coeff );
    meanMixedCoeff = mean(tempMixedCoeff, 3, 'omitnan' );
    
    tempDeformCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rDeform ).coeff );
    meanDeformCoeff = mean(tempDeformCoeff, 3, 'omitnan' );
    
    tempLocoCoeff = cat(3, postCSD_result{x}( postCSD_summary{x}.rLoco ).coeff );
    meanLocoCoeff = mean(tempLocoCoeff, 3, 'omitnan' );
    
    sp(1) = subtightplot(1,4,1,opt{:});
    if postCSD_summary{x}.nIns > 0
        plot(postCSD_opts{x}.lags,  meanInsCoeff, 'LineWidth',LW );
    end
    axis square;
    tempPos = get(gca,'Position');
    xlabel('Lag (s)'); ylabel('Coefficient'); 
    title( sprintf('Insensitive (n = %i)', postCSD_summary{x}.nIns) );
    legend(postCSD_pred{x}.name, 'Location','NorthWest', 'AutoUpdate',false)
    set(gca,'Position',tempPos);
    
    sp(2) = subtightplot(1,4,2,opt{:});
    if postCSD_summary{x}.nMixed > 0
        plot(postCSD_opts{x}.lags,  meanMixedCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Mixed (n = %i)', postCSD_summary{x}.nMixed) );
    xlabel('Lag (s)'); 
    
    sp(3) = subtightplot(1,4,3,opt{:});
    if postCSD_summary{x}.nDeform > 0
        plot(postCSD_opts{x}.lags,  meanDeformCoeff, 'LineWidth',LW ); 
    end
    axis square;
    title( sprintf('Deformation-dependent (n = %i)', postCSD_summary{x}.nDeform) );
    xlabel('Lag (s)'); %title('Deformation-dependent'); % ylabel('Coefficient');
    
    sp(4) = subtightplot(1,4,4,opt{:});
    if postCSD_summary{x}.nLoco > 0
        plot(postCSD_opts{x}.lags,  meanLocoCoeff, 'LineWidth',LW );
    end
    axis square;
    xlabel('Lag (s)'); 
    title(sprintf('Locomotion-dependent (n = %i)', postCSD_summary{x}.nLoco)); % ylabel('Coefficient'); 
    linkaxes(sp,'xy');
    
    pause;
    
    % Save the figure
    figPath = sprintf('%s%s_%s_SubtypeCoeff.tif', figDir, GLMname, expt(x).name);
    if exist(figPath,'file'), delete(figPath); end
    %print(SubtypeCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath);

    clf;
end

%% Compare pre and post-CSD results by unit
GLMname = 'preVpost_GLM_3D';

PrePostDevFig = figure('Units','normalized','OuterPosition',[0,0,1,1]);
plot( [1,2], [preCSDdevPool', postCSDdevPool'], 'color',[0,0,0,0.08] ); axis square;
hold on;
line([0.95,2.05], 0.1*[1,1], 'lineStyle','--', 'color','r');
set(gca,'Xtick',[1,2], 'XtickLabel',{'Pre','Post'});
ylabel('Deviance Explained');
title(sprintf('%i units well-fit pre-CSD,  %i post-CSD,  %i both',numel(find(preCSDdevPool > 0.1)), numel(find(postCSDdevPool > 0.1)),numel(find(preCSDdevPool > 0.1 & postCSDdevPool > 0.1))) );
figPath = sprintf('%s%s_%s_DevComparison.tif', figDir, GLMname, expt(x).name);
if exist(figPath,'file'), delete(figPath); end
print(PrePostDevFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath); %

%%
close all;
PrePostCoeffFig = figure('Units','normalized','OuterPosition',[0,0,1,1], 'color','w');
for p = 1:size(preCSDcoeffPool,2)
    subplot(1,4,p)
    JitterPlot( [preCSDcoeffPool(:,p), postCSDcoeffPool(:,p)], 0.45, 'monochrome',0.6 ); hold on;
    %plot( 1, preCSDcoeffPool(:,p), '.' ); hold on;
    %plot( 2, postCSDcoeffPool(:,p), '.' );
    plot( [1,2], [preCSDcoeffPool(:,p), postCSDcoeffPool(:,p) ], 'Color',[0,0,0,0.2] );
    set(gca,'Xtick',[1,2], 'XtickLabel',{'Pre','Post'})
    ylabel( preCSD_pred{xPresent(1)}.name{p} );
    axis square; box off;
    %xlim([0.95,2.05])
    %pause;
end
figPath = sprintf('%s%s_CoeffComparison.tif', figDir, GLMname);
if exist(figPath,'file'), delete(figPath); end
print(PrePostCoeffFig, figPath, '-dtiff', '-r300' );  fprintf('\nSaved %s\n', figPath); % 



