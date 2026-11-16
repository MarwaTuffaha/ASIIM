function [sigma,sigmaError] = ...
    sigmaCalculationTwoB(vacStep,p, ...
    yout,mwlVec,maxCoverage,population)
%sigmaCalculationFixed Caculate the vaccination rate sigma
%   returns a 4x16 matrix of sigma (one for each age and S class)

    sigmaError = p.sigmaError;

    % relabel vacStep
    i = vacStep;
    mwl = mwlVec(i);

    % --- Current cumulative second-dose coverage (fraction of population) ---
    v2Index = reshape(62*16+1:66*16,[],4);
    CV2 = (yout(end,v2Index(:,1)) ...
        + yout(end,v2Index(:,2)) ...
        + yout(end,v2Index(:,3)) ...
        + yout(end,v2Index(:,4)))./population;  % 1x16

    
    % --- Desired pop-level increase this step (kept i+1 to match your original) ---
    piDeltaCV = (maxCoverage(i+1)/100).*ones(1,16) - CV2(:)';  % 1x16 %Coverage gap
    for j = 1:16
        if piDeltaCV(j) < 0
            piDeltaCV(j) = 0;
        end
    end
    pi = piDeltaCV;  % keep naming consistent with your first-dose code

    % For second doses, only V1 can move to V2 in the ODE, so:
    v1Index = reshape(40*16+1:44*16,[],4);
    V1now = (yout(end,v1Index(:,1)) ...
           + yout(end,v1Index(:,2)) ...
           + yout(end,v1Index(:,3)) ...
           + yout(end,v1Index(:,4)));         % 1x16 counts
    whocanvaccinate = V1now;                         % 1x16 (eligible for 2nd dose)

    % Proportion to vaccinate next window scaled by eligible pool (analogous to first dose)
    tiny = 1e-12;
    piS = pi .* (population ./ max(whocanvaccinate, tiny));  % per-person prob on V1
    piS(whocanvaccinate < tiny) = 0;                         % no eligible -> no 2nd doses
    piS(piS > 0.99) = 0.99;                                  % keep finite rates
    piS(piS < 0) = 0;

    g = piS;

    % Probability -> rate over window
    sigma = zeros(4,16);
    for ii=1:16
        sigma(:,ii) = 1/mwl .* (-log(1 - g(ii))) .* ones(4,1);
    end

    % Error warning / safety
    sigma(sigma<0) = 0;
    sigma(~isfinite(sigma)) = 1/(16*7);
end
