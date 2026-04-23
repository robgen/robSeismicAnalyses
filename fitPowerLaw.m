function [a, b, sigmaNoC, powerLawWithCollapses] = fitPowerLaw(...
    IM, EDP, collapseEDP, yieldEDP)

% categorise 
collapseFLAG = zeros(size(EDP));
collapseFLAG(isinf(EDP) | EDP>=collapseEDP) = 1;

%% no-collapse regression
 
noCollapseIM    = IM(~isinf(EDP) & EDP<collapseEDP);
noCollapseEDP   = EDP(~isinf(EDP) & EDP<collapseEDP);

if yieldEDP == 0
    % single regression with all the points
    [xData, yData] = prepareCurveData( log(noCollapseIM), log(noCollapseEDP) );
    
    ft = fittype( 'poly1' );
    fitresult = fit( xData, yData, ft );
    
    a = exp(fitresult.p2);
    b = fitresult.p1;
else
    % line to interpolate the points up to yielding (get a)
    preYield = noCollapseEDP <= yieldEDP;
    
    [xData, yData] = prepareCurveData( ...
        noCollapseIM(preYield), noCollapseEDP(preYield) );
    ft = fittype( {'x'} );
    fitresult = fit( xData, yData, ft );
    a = fitresult.a;
    
    % power law of all the points, given a (get b)
    [xData, yData] = prepareCurveData( log(noCollapseIM), log(noCollapseEDP)-log(a) );
    ft = fittype( {'x'} );
    fitresult = fit( xData, yData, ft );
    b = fitresult.a;
end

resid = log(noCollapseEDP) - (log(a) + b.*log(noCollapseIM));
sigmaNoC = std( resid );
betaNoC = sigmaNoC / abs(b);

%% logistic regression for the collapse probability

IMfrag = linspace(0, max(IM), 100);
if sum(collapseFLAG) >= 3
    
    logisticParam = glmfit(log(IM), collapseFLAG, 'binomial', 'link', 'logit');
    PROBfragCollapse = 1 ./ ( 1 + exp(-logisticParam(1) -logisticParam(2)*log(IMfrag)) );
    
    powerLawWithCollapses(:,2) = ...
        (a.*IMfrag.^b) .* exp( betaNoC*(0.5./(1-PROBfragCollapse)) );
    
%     % Final approximated regression
%     forRegression = powerLawWithCollapses <= 0.10;
%     [xData, yData] = prepareCurveData( log(IMfrag(forRegression)), ...
%         log(powerLawWithCollapses(forRegression)) );
%     
%     ft = fittype( 'poly1' );
%     fitresult = fit( xData, yData, ft );
%     
%     aa = exp(fitresult.p2);
%     bb = fitresult.p1;
else
    powerLawWithCollapses(:,2) = a.*IMfrag.^b;
end
powerLawWithCollapses(:,1) = IMfrag;

end
