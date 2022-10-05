%% View CSD effects for all 3D data
ViewCSDdeform = figure('WindowState','maximized');
showVars = {'fluor', 'scaleMag', 'shearMag', 'shiftZ'};
showLims = struct('fluor',[0,5], 'scaleMag',[0,5], 'shearMag',[0,0.01], 'shiftZ',[-3,3]);

opt = {[0.01,0.01], [0.05,0.03], [0.06,0.02]};  % {[vert, horz], [bottom, top], [left, right] }
xUse = [14]; %intersect(xWave, x3D); %x3Dcsd;
spSize = [numel(xUse), numel(showVars)];
for k = 1:numel(xUse) %x3Dcsd
    zeroTick = find(csdBout{xUse(k)}(expt(xUse(k)).csd).T{1}-csdBout{xUse(k)}(expt(xUse(k)).csd).Tstart == 0);
    tempT = csdBout{xUse(k)}(expt(xUse(k)).csd).T{1}-csdBout{xUse(k)}(expt(xUse(k)).csd).Tstart;
    xTicks = 1:zeroTick-1:csdBout{xUse(k)}(expt(xUse(k)).csd).Nscan;
    for s = 1:numel(showVars)
        subtightplot(spSize(2), spSize(1),  sub2ind([numel(xUse), 4], k, s), opt{:})
        tempData = csdBout{xUse(k)}(expt(xUse(k)).csd).(showVars{s}){1};
        imagesc(tempData');
        if ~strcmpi(showVars{s}, 'fluor')
            caxis(showLims.(showVars{s}));
        else
            caxis(prctile(tempData(:), [5, 99.5]))
        end
        colormap(gca, bluewhitered);
        if k == 1
            axPos = get(gca,'Position');
            CB = colorbar('WestOutside'); CB.Label.String = showVars{s}; CB.Label.FontWeight = 'bold';
            set(gca,'Ytick',[1,round(size(tempData, 2)/2),size(tempData, 2)], 'Position',axPos, 'Xtick',xTicks);
        else
            set(gca,'Ytick',[], 'Xtick',xTicks);
        end
        if s == 1, title( sprintf('x = %i: %s', xUse(k), expt(xUse(k)).name ), 'Interpreter','none' ); end % showVars{s}
        if s == numel(showVars) 
            set(gca, 'Xtick',xTicks, 'XtickLabel',round(tempT(xTicks))); 
            xlabel('Peri-CSD time (s)');
        end
        pause;
    end
end
%impixelinfo;