function resultant = AngularResultant( Theta, Radius, varargin )
% Given a polar data set, calculate angle and magnitude, as well as normalized magnitude, of the resultant.

if nargin > 2, show = varargin{1}; else, show = false; end


[xTemp, yTemp] = pol2cart( Theta, Radius );
[resAngle, resMag] = cart2pol( sum(xTemp), sum(yTemp) );
resNorm = resMag/sum(Radius);
resultant = [resNorm, resAngle, resMag];

if show
    figure;
    polarplot( Theta, Radius, 'b.' ); hold on;
    polarplot( resAngle*[1,1], max(Radius)*[0,1], 'k' );
    title( sprintf('Normalized magnitude = %2.2f', resNorm) );
end
end

