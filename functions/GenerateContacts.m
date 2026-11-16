function [c, A, con, B] = GenerateContacts(se,simStep,kval)
%errorMsg = 'Perturbation Matrix Error';
%kval contact scaling

%debug
%sprintf('phase:%d contact%d%d%d',step,phase(1),phase(2),phase(3))

k = kval(simStep);

%% Pertubation matrices

SF = ones(16);
WF = ones(16);
OF = ones(16);

% SF = zeros(16);
% WF = zeros(16);
% OF = zeros(16);

%%
%design phase from contactModifers
phase = se.contactModifers(simStep,:);

switch phase(1)
    case 1
        SF = se.SF1;
    case 2
        SF = se.SF2;
    case 3
        SF = se.SF3;
    case 4
        SF = se.SF4;
    case 5
        SF = se.SF5;
    case 6
        SF = se.SF6;
end
        

if phase(3) == 1
    WF = se.WF1;
elseif phase(3) == 2
    WF = se.WF2;
    elseif phase(3) == 3
    WF = se.WF3;
        elseif phase(3) == 4
    WF = se.WF4;
end

if phase(2) == 1
    OF = se.OF1;
elseif phase(2) == 2
    OF = se.OF2;
elseif phase(2) == 3
    OF = se.OF3;
elseif phase(2) == 4
    OF = se.OF4;
elseif phase(2) == 5
    OF = se.OF5;
elseif phase(2) == 6
    OF = se.OF6;
elseif phase(2) == 7
    OF = se.OF7;
end

%% Population

Pop = se.population';

%%
% Scale contact matrices by population size to 
% obtain the absolute number of contacts in each group
H1 = se.H.*repmat(Pop,1,16);
S1 = se.S.*repmat(Pop,1,16);
W1 = se.W.*repmat(Pop,1,16);
O1 = se.O.*repmat(Pop,1,16);

% Multiply each contact matrix by its transpose to account for symmetries
Hmat1 = H1.*H1';
Smat1 = S1.*S1';
Wmat1 = W1.*W1';
Omat1 = O1.*O1';

% square root the result to undo the squaring effect in the previous step
% Then normalize by population size to get per capita rates
Hmat = sqrt(Hmat1)./repmat(Pop,1,16);
Smat = sqrt(Smat1)./repmat(Pop,1,16);
Wmat = sqrt(Wmat1)./repmat(Pop,1,16);
Omat = sqrt(Omat1)./repmat(Pop,1,16);


%%
% % figure(1)
% % subplot(2,2,1)
% % %imagesc(HF)
% % title(se.dateList(step,:))
% % caxis([0 2])
% % subplot(2,2,2)
% % imagesc(SF)
% % colorbar
% % title('School')
% % caxis([0 1.25])
% % subplot(2,2,3)
% % imagesc(WF)
% % colorbar
% % title('Work')
% % caxis([0 1.25])
% % subplot(2,2,4)
% % imagesc(OF)
% % colorbar
% % title('Other')
% % caxis([0 1.25])

%sprintf('SF.*Smat: %f',sum(sum(SF.*Smat)))
%sprintf('WF.*Wmat: %f',sum(sum(WF.*Wmat)))

% USc = Hmat + WF.*Wmat + OF.*Omat + SF.*Smat;

% aggregate contact matrix, accounting for contributions from all contact types
%USc = Hmat + (ones(16) - WF).*Wmat + (ones(16)-OF).*Omat + (ones(16)-SF).*Smat;
USc = Hmat + (WF).*Wmat + (OF).*Omat + ( SF).*Smat;
B = sum(USc,2); % total number of contacts for each specific group 

con = k.*USc; % contacts scaled by k

A = sum(con,2); % activity

c = con./repmat(A,1,16); % fraction of contacts

% if simStep==20
%     disp("pause")
% end



