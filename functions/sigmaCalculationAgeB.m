function [sigma,sigmaError] =...
    sigmaCalculationAgeB(vacStep,p,...
    yout,vaccage,mwlVec,maxCoverage,population)

% mwlVec is mitagationWindowLength Vector

% Caculate the vaccination rate sigma
% Returns a 4x16 matrix of sigma (one for each age and S class)
    sigmaError = p.sigmaError;

    % relabel vacStep
    i = vacStep;
    
    mwl = mwlVec(i);
 
        
    v1Index = reshape(58*16+1:62*16,16,4);
    CV1 = (yout(end,v1Index(:,1))...
    + yout(end,v1Index(:,2))...
    + yout(end,v1Index(:,3))...
    + yout(end,v1Index(:,4)))./population;

    %propPop = 1.10; %?????????????????
    propPop = 1;
    % if i+1>112
    %     disp("error")
    % end
    
    %Coverage gap
    %piDeltaCV = (maxCoverage(i+1)/100).*ones(1,16).*propPop.*p.rho - CV1(:)';
    piDeltaCV = (maxCoverage(i+1)/100).*ones(1,16).*propPop - CV1(:)';
    % if >0  more vaccinations are needed to reach the target
    % if <=0 the target is already met , and no additional vaccination is required
    for j = 1:16 % ensuring no negative vaccination demand
        if piDeltaCV(j) < 0 
            piDeltaCV(j) = 0;
        end
    end

pi = piDeltaCV;

rp = postProcess(yout(end,:));

    %Adds all individuals in Susceptible (S) and Exposed (E1, E2, E3) compartments.
    %Includes a subset of Infected (I(:,:,2)), possibly indicating a mild infection stage where vaccination is still allowed.
        % whocanvaccinate = ...
        % sum(S(:,:,1:4)+...
        % E1(:,:,1:4)+...
        % E2(:,:,1:4)+...
        % E3(:,:,1:4),3)+...
        % I(:,:,2);
        whocanvaccinate =  rp(1).S+rp(2).S+rp(3).S+rp(4).S;

    %Proportion to vaccinate next month scaled by the population of the
    %classes to vaccinate

    sigma = zeros(4,16);
    % Support both vaccine-target formats:
    %   Legacy format: a nonzero entry in column j contains the label j.
    %   Weighted format: entries contain relative weights (e.g., 0, 0.2, 1).
    % Legacy files, including the Ontario file, are converted internally to
    % binary eligibility so that their original behaviour is unchanged.
    expectedLabels = repmat(1:16,size(vaccage,1),1);
    positiveEntries = vaccage>0;
    isLegacyFormat = all(vaccage(positiveEntries) == ...
                         expectedLabels(positiveEntries));

    if isLegacyFormat
        priorityWeight = double(vaccage(i,:)>0);
    else
        priorityWeight = vaccage(i,:);
    end
    eligibleMask = priorityWeight>0;
    tiny = 1e-12;

    if any(whocanvaccinate(eligibleMask) > tiny)
        % Baseline vaccination probability used by the original code.
        piS = pi.*(population./max(whocanvaccinate,tiny));
        piS(whocanvaccinate <= tiny) = 0;

        % Redistribute, rather than increase, the intended doses. This
        % keeps the total number of doses before the 0.99 cap unchanged.
        baselineDoses = sum(piS(eligibleMask).*whocanvaccinate(eligibleMask));
        weightedPiS = piS.*priorityWeight;
        weightedDoses = sum(weightedPiS(eligibleMask).*whocanvaccinate(eligibleMask));
        if weightedDoses > tiny
            piS = weightedPiS.*(baselineDoses/weightedDoses);
        else
            piS(:) = 0;
        end

        g = piS(eligibleMask); % Filter for eligible age groups

        %%Error test - who can vaccinate is wrong
        g(g>.99)=.99; %Ensures probabilities do not exceed .99

        vaccagelist = find(eligibleMask); % Columns identify the age groups
        for ii=1:length(vaccagelist)
            sigma(:,vaccagelist(ii)) =1/mwl.*(-log(1-g(ii))).*ones(4,1);
            %-log(1 - g) is the standard transformation from probability to rate in an exponential model.
        end
    end


    %Error warrning
    
    %Blunt force error correction
    sigma(sigma<0)=0;
    
    %%Maximum vaccination rate
    maxVacRate = (sum(population)*0.0270)/16;%0.0140;
    maxSigma = sum((maxVacRate./whocanvaccinate));
    if sum(sigma(1,:))~=0
        scale = maxSigma/sum(sigma(1,:));
        if scale < 1
            sigma = sigma.*scale;
        end
    end

end
