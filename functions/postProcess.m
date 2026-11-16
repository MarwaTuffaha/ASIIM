function [r] = postProcess(yout1,N)
%POSTPROCESS Process yout form ImmunoEpi.m to a more convinent form
%r(*).ES yeilds the sum accross E^n_* for n=1..3 and age class
arguments
    yout1
    N (1,1) int8 = 16
end
N=16;

%% Susceptible Compartments

% sIndex: This index is used to organize and access the susceptible (S) compartments 
% from the model's output. It reshapes the range 1 to 8*N into an 8xN array, 
% where N represents the number of age classes. Each row in this array corresponds 
% to a different state within the susceptible population, and columns 
% represent the different age classes.

sIndex = reshape(reshape(1:8*N,N,[]),8,[]);

for i=1:8
    r(i).S = yout1(:,sIndex(i,:));
end

%% Exposed Compartments

% eIndex (from 8*N+1 to 32*N): This index range maps to the exposed compartments 
% (dE1, dE2, dE3) in the dy vector. It reshapes the specified range into an 8xN 
% array for structured access to data related to the exposed individuals in various 
% age classes, segmented into different exposure stages.

eIndex = reshape(reshape(8*N+1:32*N,N,[]),8,[]);

r(2).E1 = yout1(:,eIndex(2,1:N));
r(2).E2 = yout1(:,eIndex(3,1:N));
r(2).E3 = yout1(:,eIndex(4,1:N));
r(2).ES = sum([r(2).E1,r(2).E2,r(2).E3],2);

r(3).E1 = yout1(:,eIndex(2,N+1:2*N));
r(3).E2 = yout1(:,eIndex(3,N+1:2*N));
r(3).E3 = yout1(:,eIndex(4,N+1:2*N));
r(3).ES = sum([r(3).E1,r(3).E2,r(3).E3],2);

r(4).E1 = yout1(:,eIndex(2,2*N+1:3*N));
r(4).E2 = yout1(:,eIndex(3,2*N+1:3*N));
r(4).E3 = yout1(:,eIndex(4,2*N+1:3*N));
r(4).ES = sum([r(4).E1,r(4).E2,r(4).E3],2);

r(6).E1 = yout1(:,eIndex(6,1:N));
r(6).E2 = yout1(:,eIndex(7,1:N));
r(6).E3 = yout1(:,eIndex(8,1:N));
r(6).ES = sum([r(6).E1,r(6).E2,r(6).E3],2);

r(7).E1 = yout1(:,eIndex(6,N+1:2*N));
r(7).E2 = yout1(:,eIndex(7,N+1:2*N));
r(7).E3 = yout1(:,eIndex(8,N+1:2*N));
r(7).ES = sum([r(7).E1,r(7).E2,r(7).E3],2);

r(8).E1 = yout1(:,eIndex(6,2*N+1:3*N));
r(8).E2 = yout1(:,eIndex(7,2*N+1:3*N));
r(8).E3 = yout1(:,eIndex(8,2*N+1:3*N));
r(8).ES = sum([r(8).E1,r(8).E2,r(8).E3],2);



%% Infected Compartments

% iIndex: Infected (I) compartments from the model's output. The index range begins from
% (16+16)*N+1 to (16+16+8)*N, reflecting the position of the infected compartments in the
% dy vector after susceptible and exposed compartments. It reshapes this range into an 8xN
% array, where N is the number of age classes. Each row represents a different state within
% the infected population, and columns correspond to the different age classes. This setup is 
% used for extracting and summarizing infected data across various age groups and 
% states, with 'r(i).I' storing the data for each state and 'r(i).IS' providing 
% the sum of infected individuals in each state.

iIndex = reshape(reshape((16+16)*N+1:(16+16+8)*N,N,[]),8,[]);

for i=1:8
    r(i).I = yout1(:,iIndex(i,:));
    r(i).IS = sum(yout1(:,iIndex(i,1:N)),2);
end

%% Vaccinated Compartments

% First set of Vaccinated compartments (V1 to V4 for the first group)
% vIndex: Index range for the first group of vaccinated compartments in the dy vector.
% This range is reshaped into a 4xN array, where N represents the number of age classes.
% Each row corresponds to a different vaccination state (V1 to V4) for the first group.
vIndex = reshape(40*N+1:44*N, 4, N);
r(1).V1 = yout1(:, vIndex(1, 1:N));
r(1).V2 = yout1(:, vIndex(2, 1:N));
r(1).V3 = yout1(:, vIndex(3, 1:N));
r(1).V4 = yout1(:, vIndex(4, 1:N));

% Second set of Vaccinated compartments (V1 to V4 for the second group)
% vIndex: Adjusted to the next range in the dy vector for the second group.
vIndex = reshape(44*N+1:48*N, 4, N);
r(2).V1 = yout1(:, vIndex(1, 1:N));
r(2).V2 = yout1(:, vIndex(2, 1:N));
r(2).V3 = yout1(:, vIndex(3, 1:N));
r(2).V4 = yout1(:, vIndex(4, 1:N));

% Vaccinated with Booster (VB) for groups 1 to 4
% vIndex: Index range for vaccinated compartments with booster doses.
% This is mapped to the corresponding sections in the dy vector.
vIndex = reshape(48*N+1:52*N, 4, N);
r(1).VB = yout1(:, vIndex(1, 1:N));
r(2).VB = yout1(:, vIndex(2, 1:N));
r(3).VB = yout1(:, vIndex(3, 1:N));
r(4).VB = yout1(:, vIndex(4, 1:N));


%% Cumulative Cases and Other Compartments

% cIndex: Index range for cumulative cases (C) compartments. It reshapes the range
% from 52*N+1 to 55*N into a matrix for structured access to cumulative case data.
% Each column represents a different group.
cIndex = reshape(52*N+1:58*N, [], 6);

r(2).C = yout1(:, cIndex(:, 1));
r(3).C = yout1(:, cIndex(:, 2));
r(4).C = yout1(:, cIndex(:, 3));
r(6).C = yout1(:, cIndex(:, 4));
r(7).C = yout1(:, cIndex(:, 5));
r(8).C = yout1(:, cIndex(:, 6));

% % % %%
% % % cIndex = reshape(52*N+1:55*N,[],3);
% % % 
% % % r(2).C = yout1(:,cIndex(:,1));
% % % r(3).C = yout1(:,cIndex(:,2));
% % % r(4).C = yout1(:,cIndex(:,3));
% % % 
% % % r(1).CIS = yout1(:,end-N*14+1:end-N*13);
% % % r(2).CIS = yout1(:,end-N*13+1:end-N*12);
% % % r(3).CIS = yout1(:,end-N*12+1:end-N*11);
% % % 
% % % r(1).CIV = yout1(:,end-N*11+1:end-N*10);
% % % r(2).CIV = yout1(:,end-N*10+1:end-N*9);
% % % 
% % % r(2).red = yout1(:,end-N*8+1:end-N*7);
% % % r(3).red = yout1(:,end-N*7+1:end-N*6);
% % % r(4).red = yout1(:,end-N*6+1:end-N*5);
% % % 
% % % r(2).blue = yout1(:,end-N*5+1:end-N*4);
% % % r(3).blue = yout1(:,end-N*4+1:end-N*3);
% % % r(4).blue = yout1(:,end-N*3+1:end-N*2);
% % % 
% % % %r(2).yellow = yout1(:,end-N*3+1:end-N*2);
% % % r(3).yellow = yout1(:,end-N*2+1:end-N);
% % % r(4).yellow = yout1(:,end-N+1:end);
% % % 
v1Index = reshape(58*N + 1 : 62*N,[],4);
r(1).CV1 = yout1(:,v1Index(:,1));
r(2).CV1 = yout1(:,v1Index(:,2));
r(3).CV1 = yout1(:,v1Index(:,3));
r(4).CV1 = yout1(:,v1Index(:,4));
% % % 
% % % %r(1).CDG = yout1(:,end-N*12:end-N*11);
% % % % % 
v2Index = reshape(62*N + 1 : 66*N,[],4);
r(1).CV2 = yout1(:,v2Index(:,1));
r(2).CV2 = yout1(:,v2Index(:,2));
r(3).CV2 = yout1(:,v2Index(:,3));
r(4).CV2 = yout1(:,v2Index(:,4));

v3Index = reshape(66*N + 1 : 70*N,[],4);
r(1).CVB = yout1(:,v3Index(:,1));
r(2).CVB = yout1(:,v3Index(:,2));
r(3).CVB = yout1(:,v3Index(:,3));
r(4).CVB = yout1(:,v3Index(:,4));

% figure('Name','Cumulative Infected')
% hold on
% plot(sum(r(2).C,2))
% plot(sum(r(3).C,2))
% plot(sum(r(4).C,2))
end

