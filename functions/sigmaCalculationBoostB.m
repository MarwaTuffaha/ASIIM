function [sigma,sigmaError] =...
    sigmaCalculationBoostB(vaccineWindow,p,...
    yout,vacAge,maxCoverage,mitagationWindowLength,deltaCoverageBoost,population)
%sigmaCalculationFixed Caculate the vaccination rate sigma
%   returns a 4x16 matrix of sigma (one for each age and S class)
    sigmaError = p.sigmaError;

    % relabel monthOfVacProgram
    i = vaccineWindow;
    
    mwl = mitagationWindowLength(i);
 
        
    vBIndex = reshape(66*16+1:70*16,16,4);
    %vBIndex = reshape(52*16+1:56*16,16,4);
    CVB = (yout(end,vBIndex(:,1))...
    + yout(end,vBIndex(:,2))...
    + yout(end,vBIndex(:,3))...
    + yout(end,vBIndex(:,4)))./population;

    %propPop = 1.10;
    propPop = 1;
    
    %Coverage gap
    %piDeltaCV = (maxCoverage(i)/100).*ones(1,16).*propPop.*p.rho - CVB(:)';
    piDeltaCV = (maxCoverage(i)/100).*ones(1,16).*propPop - CVB(:)';
    
    for j = 1:16
        if piDeltaCV(j) < 0
            piDeltaCV(j) = 0;
        end
    end

pi = piDeltaCV;

rp = postProcess(yout(end,:));

whocanvaccinate = rp(5).S+rp(6).S+rp(7).S+rp(8).S;

% whocanvaccinate = rp(5).S+rp(6).S+rp(7).S+rp(8).S...
%     +rp(6).E1+rp(7).E2+rp(8).E3...
%     +rp(6).E1+rp(7).E2+rp(8).E3...
%     +rp(6).E1+rp(7).E2+rp(8).E3...
%     +rp(6).I;

% %     for ii=1:4
% %         S(:,:,ii) = monthyout5(end,1+ii-1:4:4*16);
% %         V1(:,:,ii) = monthyout5(end,40*16+1+ii-1:4:44*16);
% %         V2(:,:,ii) = monthyout5(end,44*16+1+ii-1:4:48*16);
% %         E1(:,:,ii) = monthyout5(end,8*16+1+ii-1:8:16*16);
% %         E2(:,:,ii) = monthyout5(end,16*16+1+ii-1:8:24*16);
% %         E3(:,:,ii) = monthyout5(end,24*16+1+ii-1:8:32*16);
% %         I(:,:,ii) = monthyout5(end,32*16+1+ii-1:8:40*16);
% %     end
% %     
% %         whocanvaccinate = ...
% %         sum(S(:,:,:)+...
% %         E1(:,:,:)+...
% %         E2(:,:,:)+...
% %         E3(:,:,:),3)+...
% %         I(:,:,2);

    %Proportion to vaccinate next month scaled by the population of the
    %classes to vaccinate
    % Support both vaccine-target formats:
    %   Legacy format: a nonzero entry in column j contains the label j.
    %   Weighted format: entries contain relative weights (e.g., 0, 0.2, 1).
    % Legacy files, including the Ontario file, are converted internally to
    % binary eligibility so that their original behaviour is unchanged.
    expectedLabels = repmat(1:16,size(vacAge,1),1);
    positiveEntries = vacAge>0;
    isLegacyFormat = all(vacAge(positiveEntries) == ...
                         expectedLabels(positiveEntries));

    if isLegacyFormat
        priorityWeight = double(vacAge(i,:)>0);
    else
        priorityWeight = vacAge(i,:);
    end
    eligibleMask = priorityWeight>0;
    tiny = 1e-12;
    piS = pi.*(population./max(whocanvaccinate,tiny));
    piS(whocanvaccinate <= tiny) = 0;

    % Redistribute, rather than increase, the intended doses. This keeps
    % the total number of doses before the 0.99 cap unchanged.
    baselineDoses = sum(piS(eligibleMask).*whocanvaccinate(eligibleMask));
    weightedPiS = piS.*priorityWeight;
    weightedDoses = sum(weightedPiS(eligibleMask).*whocanvaccinate(eligibleMask));
    if weightedDoses > tiny
        piS = weightedPiS.*(baselineDoses/weightedDoses);
    else
        piS(:) = 0;
    end

    g = piS(eligibleMask);
    
    %%Error test - who can vaccinate is wrong
    g(g>.99)=.99;
    
    sigma = zeros(4,16);
    vaccagelist = find(eligibleMask); % Columns identify the age groups
    for ii=1:length(vaccagelist)
        sigma(:,vaccagelist(ii)) =1/mwl.*(-log(1-g(ii))).*ones(4,1);
    end

    %Error warrning
    
    %Blunt force error correction
    sigma(sigma<0)=0;
    sigma(isnan(sigma))=0;
    
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
