function y_int = immunoEpiBoostYINT(se, population)
if ~se.useInitialConditions
infected =  se.IO;
agegroup_infected = se.ageGroupInfected;


S_int = zeros(8,se.N); %zero matrix 8x16
S_int(1,:) = se.SIC.*population(1:se.N); %SIC is 1x16 vector of ones, =population vector S_int(row one)
S_int(1,agegroup_infected) = S_int(1,agegroup_infected)-infected; %columns 5-10 contain population - infected

E1_int = zeros(8,se.N); %zero matrix 8x16
E2_int = zeros(8,se.N); %zero matrix 8x16
E3_int = zeros(8,se.N); %zero matrix 8x16

I_int = zeros(8,se.N); %zero matrix 8x16
%I_int(2,agegroup_infected) = 10/6.*infected; %I_int row 2 col.5-10 = 10/6*infected
%I_int(3,agegroup_infected) = 1/6.*infected;  %I_int row 3 col.5-10 = 1/6*infected

%I_int(2,agegroup_infected) = 0.9.*infected; %I_int row 2 col.5-10 = 0.9*infected
%I_int(3,agegroup_infected) = 0.1.*infected;  %I_int row 3 col.5-10 = 0.1*infected

for i=1:length(agegroup_infected)
    I_int(2,agegroup_infected(i)) = se.p1(1,agegroup_infected(i)).*infected; %I_int row 2 col.5-10 = 0.9*infected
    I_int(3,agegroup_infected(i)) = se.p1(2,agegroup_infected(i)).*infected;  %I_int row 3 col.5-10 = 0.1*infected
    I_int(4,agegroup_infected(i)) = se.p1(3,agegroup_infected(i)).*infected; %I_int row 2 col.5-10 = 0.9*infected
end
%% temp debug --- need zeros
V1_int = zeros(4,se.N); %zero matrix 4x16
V2_int = zeros(4,se.N); %zero matrix 4x16
VB_int = zeros(4,se.N); %zero matrix 4x16

% y_int = readmatrix('data/July15Start.xlsx', 'Range', 'A1:AVW1');
% y_int = July15Start;
y_int = [reshape(S_int,8*se.N,1);reshape(E1_int,8*se.N,1);reshape(E2_int,8*se.N,1);... %reshapes S_int,E123_int,I_int into 128x1 vector
    reshape(E3_int,8*se.N,1);reshape(I_int,8*se.N,1);...
    reshape(V1_int,4*se.N,1); reshape(V2_int,4*se.N,1); reshape(VB_int,4*se.N,1);...
    reshape(zeros(32,se.N),[],1)];
else
    y_int = se.yint';
end

end
