%% compare all/exp/comp decodings
decodingDataPath = 'D:\MATLAB\Dura\decodingAnalysis.mat';
if exist(decodingDataPath, 'file')
    fprintf('\nLoading %s', decodingDataPath); 
    load(decodingDataPath);
else
    thetaBinWidth = pi/8; % /2
    NactBin = 10;
    thetaLim = [-pi,pi; -pi,-pi/2; 0,pi/2]; %
    minFrac = 0.5; % 0.67;
    verName = '10act_50pct'; % '10act_67pct_rev';
    decodeName = {'All Angles','Expansion Only','Compression Only'};
    clearvars errStruct;
    tic
    for x = xPresent
        % Concatenate signals
        zCat = [fluor{x}.z];
        zCat = vertcat(zCat.ROI);
        scaleCat = [vertcat(deform{x}.scaleAP), vertcat(deform{x}.scaleML)];
        parfor n = 1:size(thetaLim,1)
            fprintf('\nx = %i, n = %i:', x, n );
            [~, ~, errStruct(x,n), ~, ~, ~, ~, decodeLim] = ...
                DecodeDirection(scaleCat, zCat, 'minFrac',minFrac, 'NactBin',NactBin, 'thetaLim',thetaLim(n,:), 'dTheta',thetaBinWidth, 'thetaShift',0, 'show',false); %
        end
        toc
    end
    obsMAE = reshape( [errStruct(xPresent,:).MAE], Npresent, 3 );
    nullMAEmean = reshape( [errStruct(xPresent,:).nullMAEmean], Npresent, 3 );
    nullMAEci = permute( reshape( [errStruct(xPresent,:).nullMAEci], 2, Npresent, 3 ), [2,1,3] );
    
    fprintf('\nSaving %s', decodingDataPath);
    save(decodingDataPath, 'thetaBinWidth','NactBin','thetaLim','minFrac','verName','decodeName','errStruct','obsMAE','nullMAEmean','nullMAEci');
end

FS = 13; 
MS = 14;
opt = {[0.09,0.08], [0.09,0.09], [0.09, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
close all; clearvars sp h;
MaxLikeConfInt = figure('Units','normalized', 'OuterPosition',[0, 0, 1, 1], 'Color','w'); %  
for n = 1:3
    subtightplot(3,1,n,opt{:});
    for x = xPresent
        k = find(xPresent == x);
        line( k*[1,1], errStruct(x,n).nullMAEci, 'color','k', 'LineWidth',1.5 ); hold on;
    end
    plot( 1:Npresent, [errStruct(xPresent,n).MAE], 'b.', 'MarkerSize',MS ); hold on;
    title( sprintf('%s: %i of %i experiments outperform chance', decodeName{n}, sum(obsMAE(:,n) - nullMAEci(:,1,n) < 0), Npresent) );
    set(gca,'box','off', 'tickDir','out', 'FontSize',FS, 'Xtick',1:Npresent);
    ylabel('Med. Abs. Error (deg)');
end
h(1) = plot([1,2], nan(1,2), 'color','k', 'LineWidth',1.5 );
h(2) = plot(1, NaN, 'b.', 'MarkerSize',MS ); 
legend(h, '5-95% CI', 'Decoded');
xlabel('Experiment'); 
figName = sprintf('MaxLikeComparison_ConfInt_%s.tif', verName ); % %iang_%iact  , trainData(1).bin.dir.Nbin, NactBin
figPath = ['C:\Users\ablaeser\Documents\Afferent Paper\MaxLikeResults\', figName];
print( MaxLikeConfInt, figPath, '-dtiff', '-r300' );  % 

FS = 14;
opt = {[0.05,0.08], [0.06,0.06], [0.07, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
close all;
MaxLikeErrorComp = figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w'); %  % 'WindowState','maximized'
for n = 1:3
    subtightplot(1,3,n,opt{:})
    plot( [1,2], [nullMAEci(:,1,n), obsMAE(:,n)], 'k' ); hold on;
    plot( 1, nullMAEci(:,1,n), 'k.' );
    plot( 2, obsMAE(:,n), 'k.' );
    title( sprintf('%s: %i of %i experiments outperform chance', decodeName{n}, sum(obsMAE(:,n) - nullMAEci(:,1,n) < 0), Npresent) );
    xlim([0.75, 2.25]);
    ylabel('Med. Abs. Error (deg)');
    set(gca,'box','off', 'tickDir','out', 'FontSize',FS, 'Xtick',1:2, 'XtickLabel',{'Chance','Decoded'});
    axis square;
end
figName = sprintf('MaxLikeError_Lines_%s.tif', verName ); % %iang_%iact  , trainData(1).bin.dir.Nbin, NactBin
figPath = ['C:\Users\ablaeser\Documents\Afferent Paper\MaxLikeResults\', figName];
%print( MaxLikeErrorComp, figPath, '-dtiff', '-r300' );  % 


%% Divide the data into training and testing subsets. Bin training data by angle
thetaLim = [-pi,pi]; % [0,pi/2]; %[-pi, -pi/2]; %  [0,pi/2];
thetaBinWidth = pi/6; % /2
NactBin = 10;
verName = 'exp_10act_pi6_rev1'; % 4bin_

P = repmat( struct('theta',[], 'a',[], 'a_theta',[], 'theta_a',[]), 1, Nexpt );
Pdens = repmat( struct('a',[], 'a_theta',[], 'binEdges',[]), 1, Nexpt );
maxLikeTheta = cell(1,Nexpt);  likeTheta = cell(1,Nexpt);  
clearvars checkData trainData errStruct decodeLim
for x = 1:Nexpt % par
    fprintf('\nx = %i:', x );
    [maxLikeTheta{x}, likeTheta{x}, errStruct(x), P(x), Pdens(x), trainData(x), checkData(x), decodeLim(x)] = ...
        DecodeDirection([merged(x).dScML, merged(x).dScAP], merged(x).sub, 'thetaLim',thetaLim, 'dTheta',thetaBinWidth, 'thetaShift',0, 'NactBin', NactBin, 'show',true, 'name',verName, 'pause',1); %  
end

%% Summarize the results across experiments
FS = 16; MS = 10;
opt = {[0.05,0.08], [0.06,0.06], [0.07, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
close all;
MaxLikeAccuracy = figure('Units','normalized', 'Color','w', 'WindowState','maximized'); %  %    , 'Position',[0.1, 0.05, 0.7, 0.7]
subtightplot(1,3,1, opt{:});
for x = 1:Nexpt  
    line( x*[1,1], errStruct(x).nullHRci, 'color','k' ); hold on;
    %errorbar( x, mean(errStruct(x).nullHR), std(errStruct(x).nullHR), 'ko' ); hold on; 
end
plot( 1:Nexpt, [errStruct.hitRate], 'b.', 'MarkerSize',MS ); 
H(1) = line( [1,1], errStruct(1).nullHRci, 'color','k' );
H(2) = plot( 1, errStruct(1).hitRate, 'b.', 'MarkerSize',MS);
legend(H, 'Null CI', 'Decoder');

NuseBin = sum(~cellfun( @isempty, trainData(1).bin.dir.ind ));
%line([0,Nexpt+1], (1/NuseBin)*[1,1], 'lineStyle','--', 'color','k');
ylim([0,0.3]);
axis square;
set(gca,'box','off', 'tickDir','out', 'FontSize',FS, 'Xtick',1:Nexpt);
ylabel('Hit Rate'); xlabel('Experiment');
%title( sprintf('%i angular bins, %i activity bins', trainData(1).bin.dir.Nbin, NactBin ) ) ;

subtightplot(1,3,2, opt{:});
for x = 1:Nexpt
    line( x*[1,1], errStruct(x).nullMAEci, 'color','k' ); hold on;
    %errorbar( x, mean(errStruct(x).nullMAE), std(errStruct(x).nullMAE), 'ko' ); 
end
plot( 1:Nexpt, [errStruct.MAE], 'b.', 'MarkerSize',MS ); hold on;
axis square;
set(gca,'box','off', 'tickDir','out', 'FontSize',FS, 'Xtick',1:Nexpt);
xlabel('Experiment'); ylabel('Median Absolute Error (deg)');  % ylabel('Mean Squared Error (deg^2)'); 


subtightplot(1,3,3, opt{:});
plot( [1,2], [ErrLimMix; ErrObsMix], 'color',[0,0,0,0.1] ); hold on;
h(1) = errorbar([1,2], [mean(ErrLimMix), mean(ErrObsMix)], [SEM(ErrLimMix), SEM(ErrObsMix)], 'LineWidth',2, 'color','k' );
line([1,2], 100*[1,1], 'color','k', 'LineWidth',2);
text(1.5, 105, sprintf('p = %2.3f', pMix), 'HorizontalAlignment','center' );
% Pure expansion bins
plot( [3,4], [ErrLimExp; ErrObsExp], 'color',[0,0,0,0.1] );
h(2) = errorbar([3,4], [mean(ErrLimExp), mean(ErrObsExp)], [SEM(ErrLimExp), SEM(ErrObsExp)], 'LineWidth',2, 'color','b' );
line([3,4], 100*[1,1], 'color','k', 'LineWidth',2);
text(3.5, 105, sprintf('p = %2.3f', pExp), 'HorizontalAlignment','center' );
% Pure compression
plot( [5,6], [ErrLimComp; ErrObsComp], 'color',[0,0,0,0.1] );
h(3) = errorbar([5,6], [mean(ErrLimComp), mean(ErrObsComp)], [SEM(ErrLimComp), SEM(ErrObsComp)], 'LineWidth',2, 'color','r' );
line([5,6], 100*[1,1], 'color','k', 'LineWidth',2);
xlim([0.5,6.5]); 
text(5.5, 105, sprintf('p = %2.3f', pComp), 'HorizontalAlignment','center' );
ylabel('Median Absolute Error (deg)');
set(gca,'Xtick',1:6, 'XtickLabel', {'Chance (5%)', 'Decoder', 'Chance (5%)', 'Decoder', 'Chance (5%)','Decoder'}, 'box','off');
xtickangle(45);
axis square;
legend(h, 'Mixed','Pure Expansion','Pure Compression', 'Location','best')

figName = sprintf('MaxLikeSumm_%s.tif', verName ); % %iang_%iact  , trainData(1).bin.dir.Nbin, NactBin
figPath = ['C:\Users\ablaeser\Documents\Afferent Paper\MaxLikeResults\', figName];
print( MaxLikeAccuracy, figPath, '-dtiff', '-r300' );  % 



%% Compare the accuracy of decoding for pure expansion/compression to mixed
thetaSubLim = [-pi,-pi/2; 0,pi/2];
bExp = find( checkData(1).bin.dir.cent(1:end-1) >= 0 & checkData(1).bin.dir.cent(1:end-1) <= pi/2 );
bComp = find( checkData(1).bin.dir.cent(1:end-1) >= -pi & checkData(1).bin.dir.cent(1:end-1) <= -pi/2 );
bMix = 1:checkData(1).bin.dir.Nbin; bMix([bExp,bComp]) = [];


ErrLimMix = []; ErrLimExp = []; ErrLimComp = [];
ErrObsMix = []; ErrObsExp = []; ErrObsComp = [];
for x = 1:Nexpt
    binErrLim = prctile( errStruct(x).nullMAEbin, 5, 1 );
    ErrLimMix = [ErrLimMix, binErrLim(bMix)];
    ErrLimExp = [ErrLimExp, binErrLim(bExp)];
    ErrLimComp = [ErrLimComp, binErrLim(bComp)];
    ErrObsMix = [ErrObsMix, errStruct(x).MAEbin(bMix)];
    ErrObsExp = [ErrObsExp, errStruct(x).MAEbin(bExp)];
    ErrObsComp = [ErrObsComp, errStruct(x).MAEbin(bComp)];
end

[~,pMix] = ttest( ErrLimMix, ErrObsMix, 'tail','right');
[~,pExp] = ttest( ErrLimExp, ErrObsExp, 'tail','right' );
[~,pComp] = ttest( ErrLimComp, ErrObsComp, 'tail','right' );

close all; clearvars h;
MaxLikeBinnedAccuracy = figure('color','w','Units','normalized','OuterPosition',[0.1,0.06, 0.65, 0.9]);
% Mixed bins
plot( [1,2], [ErrLimMix; ErrObsMix], 'color',[0,0,0,0.1] ); hold on;
h(1) = errorbar([1,2], [mean(ErrLimMix), mean(ErrObsMix)], [SEM(ErrLimMix), SEM(ErrObsMix)], 'LineWidth',2, 'color','k' );
line([1,2], 100*[1,1], 'color','k', 'LineWidth',2);
text(1.5, 105, sprintf('p = %2.3f', pMix), 'HorizontalAlignment','center' );
% Pure expansion bins
plot( [3,4], [ErrLimExp; ErrObsExp], 'color',[0,0,0,0.1] );
h(2) = errorbar([3,4], [mean(ErrLimExp), mean(ErrObsExp)], [SEM(ErrLimExp), SEM(ErrObsExp)], 'LineWidth',2, 'color','b' );
line([3,4], 100*[1,1], 'color','k', 'LineWidth',2);
text(3.5, 105, sprintf('p = %2.3f', pExp), 'HorizontalAlignment','center' );
% Pure compression
plot( [5,6], [ErrLimComp; ErrObsComp], 'color',[0,0,0,0.1] );
h(3) = errorbar([5,6], [mean(ErrLimComp), mean(ErrObsComp)], [SEM(ErrLimComp), SEM(ErrObsComp)], 'LineWidth',2, 'color','r' );
line([5,6], 100*[1,1], 'color','k', 'LineWidth',2);
xlim([0.5,6.5]); 
text(5.5, 105, sprintf('p = %2.3f', pComp), 'HorizontalAlignment','center' );
ylabel('Median Absolute Error (deg)');
set(gca,'Xtick',1:6, 'XtickLabel', {'Chance (5%)', 'Decoder', 'Chance (5%)', 'Decoder', 'Chance (5%)','Decoder'}, 'box','off');
xtickangle(45);
axis square;
legend(h, 'Mixed','Pure Expansion','Pure Compression', 'Location','best')
figName = 'MaxLikeBinnedAccuracy';
figPath = ['C:\Users\ablaeser\Documents\Afferent Paper\MaxLikeResults\', figName];
print( MaxLikeBinnedAccuracy, figPath, '-dtiff', '-r300' ); 


%% Summarize the results for individual experiments

figName = sprintf('MaxLikeBreakdown_%s.pdf', verName ); % %iang_%iact  , trainData(1).bin.dir.Nbin, NactBin
figPath = ['C:\Users\ablaeser\Documents\Afferent Paper\MaxLikeResults\', figName];
FS = 14;
opt = {[0.12,0.08], [0.09,0.06], [0.07, 0.04]};  % {[vert, horz], [bottom, top], [left, right] }
close all; clearvars sp;
MaxLikeBreakdown = figure('WindowState','maximized', 'color','w');
for x = find( ~isnan([errStruct.MAE]) ) % 1:Nexpt
    
    % Stretch data by direction
    dirMag = cell(thetaBin.dir.Nbin,1);
    for t = flip(2:thetaBin.dir.Nbin)
        dirMag{t,1} = testStruct(x).stretch( testStruct(x).dir >= thetaBin.dir.lim(t) & testStruct(x).dir < thetaBin.dir.lim(t+1) ); %#ok<*AGROW>
    end
    dirMag{1,1} = testStruct(x).stretch( testStruct(x).dir >= thetaBin.dir.lim(1) | testStruct(x).dir < thetaBin.dir.lim(2) ); % Account for the fact that theta.dir wraps around
    meanTestStretch = cellfun( @mean, dirMag );
    subtightplot(3,3,[1,4,7],opt{:});
    polarplot( thetaDirStretch{x}, stretchMag{x}, 'k.' ); hold on; % , 'color',[0.2,0.2,0.2,0.1]
    polarplot( decodeDir{x}, decodeMag{x}, 'r.' );
    for r = find( qualSig{x}(:,1) )', polarplot( thetaBin.dir.cent, repmat(pcResult{x}(r).thresh,thetaBin.dir.Nbin+1,1), 'color',[0,0,1,0.5] ); end
    %polarplot( thetaBin.dir.cent, repmat(testStruct(x).range(1),thetaBin.dir.Nbin+1,1), 'r', 'LineWidth',1.5 );
    %polarplot( thetaBin.dir.cent, repmat(testStruct(x).range(2),thetaBin.dir.Nbin+1,1), 'r', 'LineWidth',1.5 );
    %polarplot( thetaBin.dir.cent, [meanTestStretch; meanTestStretch(1)], 'r' ); hold off;
    rlim([0,3]);
    set(gca, 'thetatick', unique( 180 + rad2deg( trainData(x).bin.dir.lim ) ), 'FontSize', FS );
    %set(gca, 'ThetaTick',[0, 90, 180, 270], 'ThetaTickLabels',{'ML-Exp.', 'AP-Expansion','ML-Comp.', 'AP-Compression'}, 'FontSize', FS); % 
    title( sprintf('x = %i: Stretch Data, Restrictions, Binning and Thresholds', x) );
    
    sp(1) = subtightplot(3,3,2,opt{:});
    plot( likeTheta{x} );
    lgd = legend(cellstr( num2str( checkData(x).bin.dir.cent(1:end-1) ) )', 'Location',[0.65, 0.85, 0.055, 0.09]); % 
    ylabel('Likelihood of \theta'); 
    title(sprintf('Hit Rate = %2.3f  (naive = %2.3f)', errStruct(x).hitRate, errStruct(x).naiveHR));
    set(gca, 'Xtick',[], 'FontSize',FS, 'TickDir','out' )

    sp(2) = subtightplot(3,3,5,opt{:});
    plot( checkData(x).dir ); hold on;
    plot( maxLikeTheta{x} );
    legend('Data', 'Prediction', 'Location','best');
    set(gca,'Ytick',[-pi, -pi/2, 0, pi/2, pi], 'YtickLabel',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'}, 'FontSize',FS, 'Xtick',[], 'TickDir','out');
    ylim([-3.2, 3.2]);
    ylabel('\theta');

    sp(3) = subtightplot(3,3,8,opt{:});
    plot( checkData(x).mag ); hold on;
    ylabel('Stretch Magnitude'); xlabel('Test Frame');
    set(gca, 'FontSize',FS, 'TickDir','out'); % , 'Xtick',[]
    linkaxes(sp,'x');
    xlim([-Inf,Inf]);

    % Distribution of predicted and observed directions
    subtightplot(3,3,3,opt{:})
    polarhistogram( checkData(x).dir,  [-pi; checkData(x).bin.dir.lim(2:end); pi], 'normalization','probability' );
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', rad2deg(checkData(x).bin.dir.cent(1:end))+180  );
    title('Distribution of Test Stretches');
    
    subtightplot(3,3,6,opt{:});%subplot(3,3,5);
    polarhistogram( maxLikeTheta{x},  [-pi; checkData(x).bin.dir.lim(2:end); pi], 'normalization','probability' )
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', rad2deg(checkData(x).bin.dir.cent(1:end))+180  );
    title('Distribution of Predicted Stretches');
    
    % Distribution of absolute error
    subtightplot(3,3,9,opt{:})
    polarhistogram( deg2rad(errStruct(x).absErr), 6, 'Normalization','probability' );
    thetalim([0, 180]);
    title(sprintf('Absolute Error (Median = %2.2f deg)', errStruct(x).MAE));
    set(gca, 'FontSize',FS, 'TickDir','out', 'ThetaTick', rad2deg(checkData(x).bin.dir.cent(1:end))+180  );
    
    pause%(1);
    %export_fig( figPath, '-pdf', '-painters','-q101', '-append', MaxLikeBreakdown );  
    clf;
end

%% Examine probabilitiy densities
close all;
figure;
for x = 12 %:Nexpt
    for r = 1:Nclick(x)
        subplot(1,3,1);
        plot( squeeze(Ptheta{x}) );
        xlim([1,2]); ylim([0,1]);
        xlabel('Theta Bin'); ylabel('P(theta)');
        axis square;
        
        subplot(1,3,2);
        plot( Pdens(x).a{r}, 'color','k', 'lineWidth',2 ); hold on;
        for b = 1:size( Pdens(x).a_theta, 1 ),  plot( Pdens(x).a_theta{b,r} );   end
        xlabel('Activity Bin'); ylabel('P(a),  P(a|theta)');
        axis square;
        xlim([1,22]);
        title( sprintf('[x,r] = [%i, %i]', x, r ) );
        set(gca, 'FontSize', FS, 'TickDir','out', 'Xtick',[] );
        ylim([0,1]);
        
        subplot(1,3,3);
        for b = 1:size( Pdens(x).a_theta, 1 ),  plot( Pdens(x).a_theta{b,r}./Pdens(x).a{r} ); hold on; end
        xlim([1,22]);
        axis square;
        xlabel('Activity Bin'); ylabel('P(a|theta)/P(a)');
        set(gca, 'FontSize', FS, 'TickDir','out', 'Xtick',[] );
        ylim([0,2]);
        pause;
        clf;
    end
end

%% Compare data vs prediction frame by frame

figure();
for x = 1:Nexpt
    for z = 1:checkData(x).N
        polarplot( checkData(x).dir(z)*[1,1], [0,1], 'k' ); hold on;
        polarplot( maxLikeTheta{x}(z)*[1,1], [0,1], 'b' ); hold on;
        title( sprintf('z = %i:  hit = %i,  error = %2.2f', z, errStruct(x).hit(z),  errStruct(x).absErr(z) ) );  % .  , 
        legend('Data','Decoded');
        set(gca, 'thetatick', unique( 180 + rad2deg( trainData(x).bin.dir.lim ) ) );
        pause(); cla;
    end
end


%% Fix angularDist

[xFix, yFix] = pol2cart( -pi:(pi/16):pi, ones(1,numel(-pi:(pi/16):pi)) );
xFix(end+1) = 1;  yFix(end+1) = 1;
close all;
%figure; plot( xFix, yFix, '.' )
AngularDist( xFix, yFix, 'thetaRef',0, 'dThetaDir',pi/2, 'shift',pi/4, 'show',true ); %  [ fixDirData, ~, magData, fixBin, dirMag, orMag, angMagStats ] = 

 

