function AngPref = AngularPreference( Theta, Rad, OrthAng, varargin )
% Given a polar data set, quantify angular tuning using 2 methods:
%1) generate the resultant, then calculate the normalized magnitude
%2) find the maximum, along with the the orthogonal angle. Then compute index = (max-orthog)/(max+orthog)

AngPref = struct('angle',[], 'mag',[]', 'norm',[], 'pNorm',[], 'max',[], 'orth',[], 'index',[], 'pIndex',[]);
% Resultant
Resultant = AngularResultant( Theta, Rad );
AngPref.angle = Resultant.angle;
AngPref.mag = Resultant.mag ;
AngPref.norm = Resultant.norm;
%AngPref.pNorm = Resultant.pNorm;

% Get the magnitude at the angle orthogonal to the maximum
MaxOrth = AngularMaxOrth( Theta, Rad, OrthAng );
AngPref.max = MaxOrth.max;
AngPref.orth = MaxOrth.orth ;
AngPref.index = MaxOrth.index;
%AngPref.pIndex = MaxOrth.pIndex;

% Shuffle results and calculate confidence interval/p-value
if nargin > 3
    Nshuff = varargin{1};
    if Nshuff > 0
        %{
        figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        subplot(1,2,1);
        polarplot( Theta, Rad, 'k' ); hold on;
        polarplot( Resultant.angle*[1,1], Resultant.mag*[0,1], 'k' );
        polarplot( MaxOrth.max.angle, MaxOrth.max.mag, 'kx' );
        polarplot( MaxOrth.orth.angle, MaxOrth.orth.mag, 'ko' );
        title( sprintf('Normalized resultant = %1.3f   MO-Index = %1.3f', Resultant.norm, MaxOrth.index) );
        %}
        for s = flip(1:Nshuff)
            RadShuff = Rad(randperm(numel(Theta))); % Shuffle which Radius goes with which angle
            ShuffResultant(s) = AngularResultant( Theta, RadShuff ); %#ok<*AGROW>
            ShuffMaxOrth(s) = AngularMaxOrth( Theta, RadShuff, OrthAng );
            %{
            subplot(1,2,2); cla;
            polarplot( Theta, RadShuff, 'b'); hold on;
            polarplot( ShuffResultant(s).angle*[1,1], ShuffResultant(s).mag*[0,1], 'k' );
            polarplot( ShuffMaxOrth(s).max.angle, ShuffMaxOrth(s).max.mag, 'bx' );
            polarplot( ShuffMaxOrth(s).orth.angle, ShuffMaxOrth(s).orth.mag, 'bo' );
            title( sprintf('Shuffled:  Normalized resultant = %1.3f   MO-Index = %1.3f', ShuffResultant(s).norm, ShuffMaxOrth(s).index) );
            pause;
            %}
        end
        AngPref.pNorm = (numel(find([ShuffResultant.norm] > Resultant.norm)) + 1)/(Nshuff+1);
        AngPref.pIndex = (numel(find([ShuffMaxOrth.index] > MaxOrth.index)) + 1)/(Nshuff+1);
    end
end

% Show the results (optional)
if nargin < 5, show = false; else, show = varargin{2}; end
if show
    figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
    subplot(1,3,1);
    polarplot( Theta, Rad, 'b' ); hold on;
    polarplot( Resultant.angle*[1,1], Resultant.mag*[0,1], 'k' );
    polarplot( MaxOrth.max.angle, MaxOrth.max.mag, 'bx' );
    polarplot( MaxOrth.orth.angle, MaxOrth.orth.mag, 'bo' );
    title( sprintf('Normalized resultant = %1.3f (p=%1.4f)   Tuning Index = %1.3f (p=%1.4f)', Resultant.norm, Resultant.pNorm, MaxOrth.index, MaxOrth.pIndex) );
    
    subplot(1,3,2);
    ecdf([ShuffResultant.norm] ); %histogram( [ShuffResultant.norm], 'normalization', )
    hold on;
    line( Resultant.norm*[1,1], [0,1], 'color', 'k', 'LineStyle','--'); 
    xlabel('Resultant Norm. Mag.'); ylabel('CDF'); title( sprintf('Nshuff = %i', Nshuff));
    axis square;

    subplot(1,3,3);
    ecdf([ShuffMaxOrth.index] ); %histogram( [ShuffResultant.norm], 'normalization', )
    hold on;
    line( MaxOrth.index*[1,1], [0,1], 'color', 'k', 'LineStyle','--'); 
    xlabel('Max-Orth Index'); ylabel('CDF'); axis square;
    
end
end
