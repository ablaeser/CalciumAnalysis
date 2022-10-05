%% Use GLM to assess contribution of different variables
locoDeform_post_pred = cell(1,Nexpt); locoDeform_post_resp = cell(1,Nexpt); locoDeform_post_opts = cell(1,Nexpt); locoDeform_post_result = cell(1,Nexpt); locoDeform_post_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
%rLocoPreFit = cell(1,Nexpt);
GLMname = 'LocoDeform_3D_postCSD';
GLMrate = 15.49/30;
for x = x3Dcsd %xPresent %
    % GLMparallel options
    locoDeform_post_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_locoDeform_postglm');
    %locoDeform_post_opts{x}.show = true;
    locoDeform_post_opts{x}.rShow = NaN;
    locoDeform_post_opts{x}.figDir = ''; % figDir;
    locoDeform_post_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDeform_post_opts{x}.standardize = true; 
    locoDeform_post_opts{x}.trainFrac = 0.75; % 1; %
    locoDeform_post_opts{x}.Ncycle = 20;
    locoDeform_post_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDeform_post_opts{x}.CVfold = 10;
    locoDeform_post_opts{x}.nlamda = 1000;
    locoDeform_post_opts{x}.maxit = 5*10^5;
    locoDeform_post_opts{x}.minDev = 0.05;  
    locoDeform_post_opts{x}.minDevFrac = 0.1;
    locoDeform_post_opts{x}.maxP = 0.05;
    locoDeform_post_opts{x}.Nshuff = 0;  Nshuff = locoDeform_post_opts{x}.Nshuff;
    locoDeform_post_opts{x}.minShuff = 15; 
    locoDeform_post_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; % 
    locoDeform_post_opts{x}.lopo = true; %false; %
    locoDeform_post_opts{x}.frameRate = GLMrate; %expt(x).scanRate; 
    locoDeform_post_opts{x}.binSize = expt(x).scanRate/GLMrate;
    locoDeform_post_opts{x}.minShuffFrame = round( locoDeform_post_opts{x}.frameRate*locoDeform_post_opts{x}.minShuff );
    windowFrame = [ceil(locoDeform_post_opts{x}.window(1)*locoDeform_post_opts{x}.frameRate), floor(locoDeform_post_opts{x}.window(2)*locoDeform_post_opts{x}.frameRate)]; % floor(locoDeform_post_opts{x}.window*locoDeform_post_opts{x}.frameRate);
    locoDeform_post_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDeform_post_opts{x}.maxShift = max( abs(windowFrame) );
    locoDeform_post_opts{x}.Nshift = numel( locoDeform_post_opts{x}.shiftFrame );  %Nshift = locoDeform_postOpts(x).Nshift;
    locoDeform_post_opts{x}.lags = locoDeform_post_opts{x}.shiftFrame/locoDeform_post_opts{x}.frameRate;

    % Define predictors
    postRuns = expt(x).csd:expt(x).Nruns;
    tempT = BinDownMean( vertcat(Tscan{x}{postRuns}), locoDeform_post_opts{x}.binSize );
    postScans = find(tempT > csdBout{x}(expt(x).csd).Tstop + 180 & tempT < csdBout{x}(expt(x).csd).Tstop + 3600 ); % exclude the acute CSD wave from this analysis
    tempVelocityCat = BinDownMean( vertcat(loco{x}(postRuns).Vdown), locoDeform_post_opts{x}.binSize );
    tempAccelCat = BinDownMean( vertcat(loco{x}(postRuns).Adown), locoDeform_post_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    tempSpeedCat = BinDownMean( vertcat(loco{x}(postRuns).speedDown), locoDeform_post_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(postRuns).stateDown), locoDeform_post_opts{x}.binSize );
    locoDeform_post_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % LocoDeformPred(x).data = []; %LocoDeformPred(x).data = LocoDeformCat.T; LocoDeformPred(x).name = {'Time'};
    locoDeform_post_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % tempSpeedCat,, 
    locoDeform_post_pred{x}.data = locoDeform_post_pred{x}.data(postScans,:);
    locoDeform_post_pred{x}.name = {'Velocity','Acceleration','State'}; % ,'Speed',
    locoDeform_post_pred{x}.N = size(locoDeform_post_pred{x}.data, 2);
    for p = flip(1:locoDeform_post_pred{x}.N), locoDeform_post_pred{x}.lopo.name{p} = ['No ',locoDeform_post_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    locoDeform_post_pred{x}.fam.col = {1:2, 3}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDeform_post_pred{x}.fam.N = numel(locoDeform_post_pred{x}.fam.col); 
    locoDeform_post_pred{x}.fam.name = {'Kinematics', 'State'}; %{'All'};%  'Onset Time',
    
    % Define responses
    tempTrans = vertcat(deform{x}(postRuns).transMag);
    tempTrans = BinDownMean( mean( tempTrans(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_post_opts{x}.binSize ); % (:,segParams{x}.zProj)
    tempTransSpd = vertcat(deform{x}(postRuns).DtransMag);
    tempTransSpd = BinDownMean( mean( tempTransSpd(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_post_opts{x}.binSize ); % vertcat(deform{x}(postRuns).DtransMag)
    tempScaleMag = vertcat(deform{x}(postRuns).scaleMag);
    tempScaleMag = BinDownMean( mean( tempScaleMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_post_opts{x}.binSize );
    tempStretchMag = vertcat(deform{x}(postRuns).stretchMag);
    tempStretchMag = BinDownMean( mean( tempStretchMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_post_opts{x}.binSize );
    tempShearMag = vertcat(deform{x}(postRuns).shearMag);
    tempShearMag = BinDownMean( mean( tempShearMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_post_opts{x}.binSize );
    tempShearRate = vertcat(deform{x}(postRuns).DshearMag);
    tempShearRate = BinDownMean( mean( tempShearRate(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_post_opts{x}.binSize );
    tempShift = vertcat(deform{x}(postRuns).shiftZ);
    tempShiftMean =  BinDownMean( mean( tempShift(:,4:end-3), 2, 'omitnan' ), locoDeform_post_opts{x}.binSize );

    locoDeform_post_resp{x}.data = [tempTrans, tempScaleMag, tempShearMag, tempTransSpd, tempStretchMag, tempShearRate, tempShiftMean]; % tempScaleMag; %  , tempDshearMag tempExp, tempComp,  , tempCompressMean
    locoDeform_post_resp{x}.data = locoDeform_post_resp{x}.data(postScans,:);
    tempNormData = normalize(locoDeform_post_resp{x}.data, 1); % locoDeform_post_resp{x}.data
    locoDeform_post_resp{x}.data(abs(tempNormData) > 5) = NaN; % suppress extreme outliers   locoDeform_post_resp{x}.data
    locoDeform_post_resp{x}.name = {'|Translation|', '|Scale|', '|Shear|', 'Trans. Speed', '|Stretch|', 'Shear Rate', 'Z Shift'}; % , 'Z Shift', 'Compression' , 'Z Compression' 
    locoDeform_post_resp{x}.N = size(locoDeform_post_resp{x}.data, 2); 
    
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDeform_post_pred{x}.data, locoDeform_post_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDeform_post_pred{x}.data(nanFrame,:) = []; locoDeform_post_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDeform_post_opts{x}.load = true; % false; % 
    locoDeform_post_opts{x}.saveRoot = expt(x).dir; % sprintf('%s%s', expt(x).dir, locoDeform_post_opts{x}.name  ); % ''; %  
    [locoDeform_post_result{x}, locoDeform_post_summary{x}, locoDeform_post_opts{x}, locoDeform_post_pred{x}, locoDeform_post_resp{x}] = GLMparallel(locoDeform_post_pred{x}, locoDeform_post_resp{x}, locoDeform_post_opts{x}); % LocoDeform_result{x}, LocoDeform_summary(x), ~, LocoDeformPred(x), LocoDeformResp(x)
    %ViewGLM(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_post_opts{x}, LocoDeform_result{x}, LocoDeform_summary{x}); %GLMresultFig = 
    %pause;
end

%% View GLM results
for x = x3Dcsd
    locoDeform_post_opts{x}.rShow = 1:locoDeform_post_resp{x}.N; % 2; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDeform_post_opts{x}.xVar = 'Time';
    ViewGLM(locoDeform_post_pred{x}, locoDeform_post_resp{x}, locoDeform_post_opts{x}, locoDeform_post_result{x}, locoDeform_post_summary{x}); %GLMresultFig = 
end

%%
devSummMat = cell(1,Nexpt);
GLMresultFig = figure('WindowState','maximized', 'color','w');
for x = x3Dcsd
    cla;
    devSummMat{x} = [locoDeform_post_summary{x}.dev; locoDeform_post_summary{x}.lopo.dev; locoDeform_post_summary{x}.lofo.dev];
    % {
    imagesc( devSummMat{x}' );
    set(gca,'XtickLabel',  [{'All'}; locoDeform_post_summary{xPresent(1)}.lopo.name(:); locoDeform_post_summary{xPresent(1)}.lofo.name(:)], ...
        'Ytick',1:locoDeform_post_resp{xPresent(1)}.N , 'YtickLabel',locoDeform_post_resp{xPresent(1)}.name, 'TickDir','out', 'TickLength',[0.003,0]);
    title(sprintf('%s: GLM Fit Summary', locoDeform_post_opts{x}.name), 'Interpreter','none');
    xtickangle(30);
    axis image;
    CB = colorbar; CB.Label.String = 'Deviance Explained';
    caxis([0,0.2]);
    impixelinfo;
    pause;
    %}
end

devSummCat = cat(3, devSummMat{:});
devSummMed = median(devSummCat, 3, 'omitnan' );

locoDeform_post_GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; locoDeform_post_summary{x}.lopo.name(:); locoDeform_post_summary{x}.lofo.name(:)], ...
    'YtickLabel',locoDeform_post_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_post_GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Show relative explanatory value from submodels
exptColor = distinguishable_colors(Npresent);
locoDeform_post_relExpFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_post_resp{x3Dcsd(1)}.N;  
%tiledlayout(1, Ncol);
opt = {[0.05,0.04], [0.06,0.02], [0.04,0.02]};
Ndrops = locoDeform_post_pred{x3Dcsd(1)}.N+locoDeform_post_pred{x3Dcsd(1)}.fam.N;
for col = 1:Ncol
    tempDevFracMat = nan( Npresent, Ndrops );
    k = 0;
    for x = x3Dcsd
        k = k + 1;
        if locoDeform_post_result{x}(col).dev > locoDeform_post_opts{x}.minDev
            tempDevFracMat(k,:) = [locoDeform_post_result{x}(col).lopo.devFrac, locoDeform_post_result{x}(col).lofo.devFrac];
        end
    end
    subtightplot(1,Ncol,col,opt{:}); %nexttile
    JitterPlot( 1-tempDevFracMat, 0.45 ); hold on;
    line(0.5+[0,Ndrops], [0,0], 'color','k'); 
    set(gca, 'Xtick',1:Ndrops); %, 'Xticklabel'
    ylim([-0.4, 1.1]); xlim(0.5+[0,Ndrops])
    %colororder(exptColor)
    axis square;
    if col == 1, ylabel('Relative Explanatory Value'); end
    title( locoDeform_post_resp{x3Dcsd(1)}.name{col} );
    set(gca, 'XtickLabel', [locoDeform_post_pred{x3Dcsd(1)}.name, locoDeform_post_pred{x3Dcsd(1)}.fam.name])  % [locoDeform_post_summary{x}.lopo.name(:); locoDeform_post_summary{x}.lofo.name(:)]
    xtickangle(45)
end
figPath = sprintf('%s%s_ExpValSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_post_relExpFig, figPath, '-dtiff');

%% Show time filters per response
exptColor = distinguishable_colors(Npresent);
locoDeform_post_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_post_resp{x3Dcsd(1)}.N;   Nrow = locoDeform_post_pred{x3Dcsd(1)}.N;
tiledlayout(Nrow, Ncol);
k = 0;
for row = 1:Nrow
    for col = 1:Ncol
        nexttile
        tempCoeffMat = zeros( numel(locoDeform_post_opts{x3Dcsd(1)}.lags), Npresent );
        k = 0;
        for x = x3Dcsd
            k = k + 1;
            if locoDeform_post_result{x}(col).dev > locoDeform_post_opts{x}.minDev
                tempCoeffMat(:,k) = locoDeform_post_result{x}(col).coeff(:,row);
            end
        end
        plot( locoDeform_post_opts{x}.lags, tempCoeffMat ); hold on;
        colororder(exptColor)
        axis square;
        if col == 1, ylabel(sprintf('%s Coefficient', locoDeform_post_pred{x3Dcsd(1)}.name{row} )); end
        if row == 1, title( locoDeform_post_resp{x3Dcsd(1)}.name{col} ); end
        if row == Nrow, xlabel('Response Lag (s)'); end
        %pause;
    end
end

figPath = sprintf('%s%s_CoeffSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_post_GLMcoeffFig, figPath, '-dtiff');



%% Compare coefficients for same response pre/post CSD
figure('WindowState','maximized');
opt = {[0.07,0.07], [0.2,0.07], [0.1,0.04]};  % {[vert, horz], [bottom, top], [left, right] }
for x = x3Dpost
    for r = unique([locoDeform_pre_summary{x}.rGood, locoDeform_post_summary{x}.rGood]) %1:locoDeform_post_resp{x}.N 
        for p = flip(1:locoDeform_post_pred{x}.N)
            subtightplot(1,locoDeform_post_pred{x}.N, p, opt{:}); cla;
            line(locoDeform_post_opts{x}.window, [0,0], 'color','k');
            hold on;
            plot( locoDeform_post_opts{x}.lags', [locoDeform_pre_result{x}(r).coeff(:,p), locoDeform_post_result{x}(r).coeff(:,p)]); 
            xlim(locoDeform_post_opts{x}.window);
            axis square;
            ylabel(sprintf('%s Coefficient', locoDeform_pre_pred{x}.name{p}));
        end
        title( locoDeform_pre_resp{x}.name{r} );
        pause;
    end
end


