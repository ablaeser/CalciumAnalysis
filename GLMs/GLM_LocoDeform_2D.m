%% Use GLM to assess contribution of different variables
locoDeform_pred = cell(1,Nexpt); locoDeform_resp = cell(1,Nexpt); locoDeform_opts = cell(1,Nexpt); locoDeform_result = cell(1,Nexpt); locoDeform_summary = cell(1,Nexpt);
glmFigDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( glmFigDir )
GLMname = 'LocoDeform_2D_60s_accelMag'; % _noState
GLMrate = 15.49/3;
for x = xPresent %x3Dcsd
    % GLMparallel options
    locoDeform_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_locoDeformglm');
    %locoDeform_opts{x}.show = true;
    locoDeform_opts{x}.rShow = NaN;
    locoDeform_opts{x}.glmFigDir = ''; % glmFigDir;
    locoDeform_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDeform_opts{x}.standardize = true; 
    locoDeform_opts{x}.trainFrac = 0.75; % 1; %
    locoDeform_opts{x}.Ncycle = 20;
    locoDeform_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDeform_opts{x}.CVfold = 10;
    LocoDeform.opts(x).nlamda = 1000;
    LocoDeform.opts(x).maxit = 5*10^5;
    locoDeform_opts{x}.minDev = 0.05;  minDev = locoDeform_opts{x}.minDev;
    locoDeform_opts{x}.minDevFrac = 0.1;
    locoDeform_opts{x}.maxP = 0.05;
    locoDeform_opts{x}.Nshuff = 0;  Nshuff = locoDeform_opts{x}.Nshuff;
    locoDeform_opts{x}.minShuff = 15; 
    locoDeform_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; % 
    locoDeform_opts{x}.lopo = true; %false; %
    locoDeform_opts{x}.frameRate = GLMrate; %expt(x).scanRate; 
    locoDeform_opts{x}.binSize = expt(x).scanRate/GLMrate;
    locoDeform_opts{x}.minShuffFrame = round( locoDeform_opts{x}.frameRate*locoDeform_opts{x}.minShuff );
    windowFrame = [ceil(locoDeform_opts{x}.window(1)*locoDeform_opts{x}.frameRate), floor(locoDeform_opts{x}.window(2)*locoDeform_opts{x}.frameRate)]; % floor(locoDeform_opts{x}.window*locoDeform_opts{x}.frameRate);
    locoDeform_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDeform_opts{x}.maxShift = max( abs(windowFrame) );
    locoDeform_opts{x}.Nshift = numel( locoDeform_opts{x}.shiftFrame );  %Nshift = locoDeformOpts(x).Nshift;
    locoDeform_opts{x}.lags = locoDeform_opts{x}.shiftFrame/locoDeform_opts{x}.frameRate;

    % Concatenate and bin input variables 
    % Define predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).Vdown), locoDeform_opts{x}.binSize );
    tempAccelCat = BinDownMean( abs(vertcat(loco{x}(expt(x).preRuns).Adown)), locoDeform_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    tempStateCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).stateDown), locoDeform_opts{x}.binSize );
    %tempSpeedCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).speedDown), locoDeform_opts{x}.binSize );
    locoDeform_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); 
    locoDeform_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % % tempSpeedCat,,  , tempStateCat
    locoDeform_pred{x}.name = {'Velocity','Accel','State'}; % {'State'}; %,'Speed', 
    locoDeform_pred{x}.N = size(locoDeform_pred{x}.data,2);
    for p = flip(1:locoDeform_pred{x}.N), locoDeform_pred{x}.lopo.name{p} = ['No ',locoDeform_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    locoDeform_pred{x}.fam.col = {1:2, 3}; % {}; %  %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDeform_pred{x}.fam.N = numel(locoDeform_pred{x}.fam.col); 
    locoDeform_pred{x}.fam.name = {'Kinematics', 'State'}; %{};% {}; %{'All'};%  'Onset Time',

    % Define responses
    tempTrans = BinDownMean( mean( vertcat(deform{x}(expt(x).preRuns).transMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempTransSpeed = BinDownMean( mean( vertcat(deform{x}(expt(x).preRuns).DtransMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempScaleMag = BinDownMean( mean( vertcat(deform{x}(expt(x).preRuns).scaleMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempStretchMag = BinDownMean( mean( vertcat(deform{x}(expt(x).preRuns).stretchMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempShearMag = BinDownMean( mean( vertcat(deform{x}(expt(x).preRuns).shearMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempDshearMag = BinDownMean( mean( vertcat(deform{x}(expt(x).preRuns).DshearMag), 2, 'omitnan'), locoDeform_opts{x}.binSize ); %circshift(tempAccelCat, 5); %   
    locoDeform_resp{x}.data = [tempTrans, tempTransSpeed, tempScaleMag, tempStretchMag, tempShearMag, tempDshearMag]; % tempScaleMag; % , tempExp, tempComp  , tempMixed
    locoDeform_resp{x}.data = normalize(locoDeform_resp{x}.data, 1);
    locoDeform_resp{x}.data(abs(locoDeform_resp{x}.data) > 5) = NaN; % suppress extreme outliers
    locoDeform_resp{x}.name = {'|Trans|','Trans Speed','|Scale|','|Stretch|','|Shear|','|Shear| Rate'}; % ,'Exp Str','Comp Str', 'Z Shift', 'Compression'  ,'Mix Str'
    locoDeform_resp{x}.N = size(locoDeform_resp{x}.data,2); 
    
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDeform_pred{x}.data, locoDeform_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDeform_pred{x}.data(nanFrame,:) = []; locoDeform_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDeform_opts{x}.load = true; %  false; %  
    locoDeform_opts{x}.saveRoot = expt(x).dir; %sprintf('%s%s', expt(x).dir  ); %''; % , locoDeform_opts{x}.name

    [locoDeform_result{x}, locoDeform_summary{x}, locoDeform_opts{x}, locoDeform_pred{x}, locoDeform_resp{x}] = GLMparallel(locoDeform_pred{x}, locoDeform_resp{x}, locoDeform_opts{x}); 
    %ViewGLM(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_opts{x}, LocoDeform_result{x}, LocoDeform_summary{x}); %GLMresultFig = 
    %pause;
end

%% View GLM results
for x = xPresent
    locoDeform_opts{x}.figDir = '';
    locoDeform_opts{x}.rShow = 3; %locoDeform_summary{x}.rGood; %1:7; %1:LocoDeform_resp{x}.N; %NaN; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDeform_opts{x}.xVar = 'Time';
    ViewGLM(locoDeform_pred{x}, locoDeform_resp{x}, locoDeform_opts{x}, locoDeform_result{x}, locoDeform_summary{x}); %GLMresultFig = 
end

%round(30*expt(x).scanRate)
%movParam.scalebar = MakeScaleBar( [round(30*expt(x).scanRate),0], {[0,expt(X).Ncol-segParams(X).edges(1)-segParams(X).edges(2)]+0.5, [0,expt(X).Nrow-segParams(X).edges(3)-segParams(X).edges(4)]+0.5},...
%    [0.1,0.95], [0,0], 'label',false, 'color','w', 'show',false );
%%
devSummMat = cell(1,Nexpt);
%GLMresultFig = figure('WindowState','maximized', 'color','w');
for x = xPresent
    
    devSummMat{x} = [locoDeform_summary{x}.dev; locoDeform_summary{x}.lopo.dev; locoDeform_summary{x}.lofo.dev];
    %{
    cla;
    imagesc( devSummMat{x}' );
    set(gca,'XtickLabel',  [{'All'}; locoDeform_summary{xPresent(1)}.lopo.name(:); locoDeform_summary{xPresent(1)}.lofo.name(:)], ...
        'Ytick',1:locoDeform_resp{xPresent(1)}.N , 'YtickLabel',locoDeform_resp{xPresent(1)}.name, 'TickDir','out', 'TickLength',[0.003,0]);
    title(sprintf('%s: GLM Fit Summary', locoDeform_opts{x}.name), 'Interpreter','none');
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

locoDeform_GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; locoDeform_summary{x}.lopo.name(:); locoDeform_summary{x}.lofo.name(:)], ...
    'YtickLabel',locoDeform_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 2D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', glmFigDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_GLMdevFig, figPath, '-dtiff');
impixelinfo;


locoDeform_GLMdevFig = figure('WindowState','maximized', 'color','w');
tempSumm = [locoDeform_summary{xPresent}];
tempDev = vertcat(tempSumm.dev);
JitterPlot(tempDev, 0.5 ); hold on;
line([0, locoDeform_resp{x}.N]+0.5, locoDeform_opts{x}.minDev*[1,1], 'lineStyle','--', 'color','r');
axis square;
ylim([0,0.4]);
set(gca, 'Xtick',1:locoDeform_resp{x}.N, 'XtickLabel',locoDeform_resp{x}.name );
xtickangle(45);
ylabel('Deviance explained');

%% Show relative explanatory value from submodels
exptColor = distinguishable_colors(Npresent);
locoDeform_relExpFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_resp{xPresent(1)}.N;  
%tiledlayout(1, Ncol);
opt = {[0.05,0.04], [0.06,0.02], [0.04,0.02]};
Ndrops = locoDeform_pred{xPresent(1)}.N+locoDeform_pred{xPresent(1)}.fam.N;
for col = 1:Ncol
    tempDevFracMat = nan( Npresent, Ndrops );
    k = 0;
    for x = xPresent
        k = k + 1;
        if locoDeform_result{x}(col).dev > locoDeform_opts{x}.minDev
            tempDevFracMat(k,:) = [locoDeform_result{x}(col).lopo.devFrac, locoDeform_result{x}(col).lofo.devFrac];
        end
    end
    subtightplot(1,Ncol,col,opt{:}); %nexttile
    JitterPlot( 1-tempDevFracMat, 0.45 ); hold on;
    line(0.5+[0,Ndrops], [0,0], 'color','k'); 
    set(gca, 'Xtick',1:Ndrops); %, 'Xticklabel'
    ylim([-0.2, 0.8]); xlim(0.5+[0,Ndrops])
    %colororder(exptColor)
    axis square;
    if col == 1, ylabel('Relative Explanatory Value'); end
    title( locoDeform_resp{xPresent(1)}.name{col} );
    set(gca, 'XtickLabel', [locoDeform_pred{xPresent(1)}.name, locoDeform_pred{xPresent(1)}.fam.name])  % [locoDeform_summary{x}.lopo.name(:); locoDeform_summary{x}.lofo.name(:)]
    xtickangle(45)
end
figPath = sprintf('%s%s_ExpValSummary.tif', glmFigDir, GLMname );
%fprintf('\nSaving %s\n', figPath);
%print( locoDeform_relExpFig, figPath, '-dtiff');

%% Show time filters per response
exptColor = distinguishable_colors(Npresent);
locoDeform_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_resp{xPresent(1)}.N;   Nrow = locoDeform_pred{xPresent(1)}.N;
tiledlayout(Nrow, Ncol);
k = 0;
for row = 1:Nrow
    for col = 1:Ncol
        nexttile
        tempCoeffMat = nan( numel(locoDeform_opts{xPresent(1)}.lags), Npresent );
        k = 0;
        for x = xPresent
            k = k + 1;
            if locoDeform_result{x}(col).dev > locoDeform_opts{x}.minDev
                tempCoeffMat(:,k) = locoDeform_result{x}(col).coeff(:,row);
            end
        end
        plot( locoDeform_opts{x}.lags, tempCoeffMat ); hold on;
        colororder(exptColor)
        axis square;
        if col == 1, ylabel(sprintf('%s Coefficient', locoDeform_pred{xPresent(1)}.name{row} )); end
        if row == 1, title( locoDeform_resp{xPresent(1)}.name{col} ); end
        if row == Nrow, xlabel('Response Lag (s)'); end
        %pause;
    end
end

figPath = sprintf('%s%s_CoeffSummary.tif', glmFigDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_GLMcoeffFig, figPath, '-dtiff');

%% Which experiment has the best overall fit?

V = [1,3,5]
for x = xPresent
    fprintf('\nx = %i: sum = %2.2f, product = %2.4f', x, sum(locoDeform_summary{x}.dev(V)), prod(locoDeform_summary{x}.dev(V)))
end

