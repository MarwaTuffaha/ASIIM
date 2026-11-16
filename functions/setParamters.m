function [p] = setParamters(obj,sigma)
% This function sets various epidemiological and demographic parameters 
% for a COVID-19 simulation model. Parameters are stored in a structure
% 'p' for faster data access and better performance.

p.N = obj.N;

p.alpha = obj.alpha;	% susceptibility
p.beta = obj.beta;      % infectivity
p.omega = obj.omegaElement;	% waning
p.delta = obj.delta;	% disease-induced mortality (=0)
p.gamma = obj.gamma;	% recovery
p.A = obj.A;            % contacts

% Vaccination-related parameters (initial rates only, changed in sim)
p.sigma = sigma(:,:,1);    % vaccination rate - first dose
p.sigma2 = sigma(:,:,2);   % vaccination rate - second dose
p.sigmaB = sigma(:,:,3);   % vaccination rate - booster dose

% Natural events
p.mu = obj.mu;   % Natural death rate
p.tau = obj.tau; % Aging rate within compartments

% Mixing and other parameters    
p.c = obj.c;            % mixing

p.p1 = obj.p1;          % infection fraction to I2 (size(3,16)) Infection only occurs from S1, S2 and S3
p.p2 = obj.p2;          % infection fraction to I3
p.p3 = obj.p3;     

p.qq = obj.qq; % Vaccine protection from severe disease

p.qq_SLabel=obj.qq_SLabel;

p.rho = obj.rho; % Vaccine efficacy
p.epsilon = obj.epsilon; % vaccine protection from single dose

p.f = obj.f;
p.kappa = obj.kappa;
p.sigmaError = obj.sigmaError;

% Threshold for setting compartment values to zero (for numerical stability)
p.threshold = obj.threshold;


% susceptible-saturation factor label (Marwa) 0=off, 1=on
p.S_fLabel=obj.S_fLabel;
end