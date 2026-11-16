function [] = verticalLines(se,dateNS,vertAlign,fontSize,lineWidth,labels,vac)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    arguments
        se
        dateNS
        vertAlign = 'bottom';
        fontSize = 15;
        lineWidth = 2.5;
        labels = true;
        vac = true;
    end

    vacDates=se.dateNumList(se.vaccinate);
    coverageStartInd=find(se.vaccineCoverage(:,1),1,'first');
    coverageStartDate=vacDates(coverageStartInd);
% predictionstart = dateNS(1031); 
endData = dateNS(end);
%vacStart = max(se.vaccinationStart,coverageStartDate);
vacStart = se.vaccinationStart;
% caseDataEnd = dateNS(137); %For CA

if sum(se.vaccinate)==0
    vac=0;
end

if vac
        xl(1) = xline(vacStart,'--','color',[.2 .2 .2]);
        xl(1).LineWidth = lineWidth;
        xl(1).FontSize = fontSize-2;
        xl(1).LabelHorizontalAlignment = 'left';
        xl(1).LabelVerticalAlignment = vertAlign;
        if se.cut_Omicron~=1 && ~isempty(se.vaccinationBoostStart) && se.vaccinationBoostStart~=se.dateNumList(1)
        xl(4) = xline(se.vaccinationBoostStart,'--','color',[.2 .2 .2]);
        xl(4).LineWidth = lineWidth;
        xl(4).FontSize = fontSize-2;
        xl(4).LabelHorizontalAlignment = 'right';
        xl(4).LabelVerticalAlignment = vertAlign;
        xl(4).Label = 'booster dose start';
        end
end
        
%         xl(2) = xline(endData,'k:');
%         xl(2).FontSize = fontSize-2;
%         xl(2).LineWidth = lineWidth;
%         xl(2).LabelHorizontalAlignment = 'left';
%         xl(2).LabelVerticalAlignment = vertAlign;
        
        %xl(3) = xline(endData,'k:');
%         xl(3) = xline(predictionstart,'k:'); %CA
%         xl(3).FontSize = fontSize-2;
%         xl(3).LineWidth = lineWidth;
%         xl(3).LabelHorizontalAlignment = 'left'; %was right
%         xl(3).LabelVerticalAlignment = vertAlign;

        % xl(5) = xline(endData,'k:');
        % xl(5).FontSize = fontSize-2;
        % xl(5).LineWidth = lineWidth;
        % xl(5).LabelHorizontalAlignment = 'left';
        % xl(5).LabelVerticalAlignment = 'middle';
        
%         xl(5) = xline(caseDataEnd,'k:'); %For CA
%         xl(5).FontSize = fontSize-2;
%         xl(5).LineWidth = lineWidth;
%         xl(5).LabelHorizontalAlignment = 'left';
%         xl(5).LabelVerticalAlignment = 'top';
        
        if labels
            if vac
            xl(1).Label = 'vacc. start';
            if se.cut_Omicron~=1&& ~isempty(se.vaccinationBoostStart) && se.vaccinationBoostStart~=se.dateNumList(1)
                xl(4).Label = 'boost. start';
            end
            end
%            xl(2).Label = 'past';
%            xl(3).Label = 'prediction';
%            xl(5).Label = datestr(endData);
        end

end

