function [ resultant, CI ] = AngularActivity( stretchAngle, stretchMag, subActivity, varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


if nargin > 3 
    TestRange = varargin{1}; 
    RangeInd = find( stretchMag > TestRange(1) & stretchMag < TestRange(2) ); % indices of stretches within the test range
    stretchAngleRange = stretchAngle( RangeInd );
    ActRange = subActivity( RangeInd );
else
    stretchAngleRange = stretchAngle;
    ActRange = subActivity;
end

resVec = AngularResultant( stretchAngleRange, ActRange, true );
resultant.norm = resVec(1); 
resultant.angle = resVec(2); 
resultant.mag = resVec(3);

% Bootstrap confidence intervals
bootConf = bootci( 1000, @AngularResultant, stretchAngleRange, ActRange );
CI.norm = [bootConf(1,1), bootConf(2,1)];
CI.angle = [bootConf(1,2), bootConf(2,2)];
CI.mag = [bootConf(1,3), bootConf(2,3)];

opt = {[0.10,0.05], [0.10,0.1], [0.03,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
subtightplot(1,2,1,opt{:})
polarplot( stretchAngle, subActivity, '.' ); title('All Data');

subtightplot(1,2,2,opt{:})
polarplot( stretchAngleRange, ActRange, '.' ); hold on;
polarplot( resultant.angle*[1,1], resultant.mag*[0,1] );
polarplot( CI.angle([1,2]), CI.mag([1,1]), 'k--' )
polarplot( CI.angle([1,2]), CI.mag([2,2]), 'k--' )
polarplot( CI.angle([1,1]), CI.mag([1,2]), 'k--' )
polarplot( CI.angle([2,2]), CI.mag([1,2]), 'k--' )
title( sprintf('Test-Range Data Only. Normalized resultant = %2.3f, conf interval = [%1.3f, %1.3f]', resultant.norm, CI.norm(1), CI.norm(2) ) );


end

