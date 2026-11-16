function [time, prevalence] = getPrevLognormal(obj, q1LOS, medianLOS, q3LOS, K)
% getPrevLognormal
% -------------------------------------------------------------------------
% Estimate daily hospital prevalence from daily hospital admissions using a
% lognormal length-of-stay (LOS) distribution.
%
% IDEA:
%   If A(t) is the number of new admissions on day t, and p(k) is the
%   probability that a patient is still hospitalized k days after
%   admission, then
%
%       prevalence(t) = sum_{k=0}^{K-1} A(t-k) * p(k)
%
%   In words:
%     - look back k days,
%     - take the admissions on that earlier day,
%     - multiply by the fraction expected to still be in hospital,
%     - then sum over all earlier admission cohorts.
%
% INPUTS:
%   obj        - model object used by getPrediction(obj)
%   medianLOS  - empirical median LOS
%   q1LOS      - empirical 25th percentile LOS
%   q3LOS      - empirical 75th percentile LOS
%   K          - maximum number of days after admission to track
%
% OUTPUTS:
%   time       - time vector returned by getPrediction(obj)
%   prevalence - estimated daily hospital prevalence
%
% EXAMPLE:
%   % Ontario first-wave example:
%   [time, prev] = getPrevLognormal(obj, 3.6, 8.7, 18.9, 60);
%
% IMPORTANT NOTE ABOUT DISCRETIZATION:
%   LOS is modeled here as a CONTINUOUS random variable, but admissions and
%   prevalence are available as DAILY counts. Therefore, we must discretize
%   the survival function.
%
%   In this code, we use:
%
%       p(k) = P(LOS >= k)
%
%   for integer k = 0,1,2,...
%
%   That is a reasonable and commonly used choice, but it is not the only
%   possible convention. A slightly different discretization could be used
%   depending on how the hospital prevalence data are defined:
%
%     - start-of-day census,
%     - end-of-day census,
%     - daily average number of occupied beds.
%
%   So the overall method is correct, but the exact expression used for p(k)
%   is tied to a timing convention for daily prevalence.
%
%   For example, one could alternatively use:
%
%       p(k) = P(LOS >= k + 0.5)
%
%   if a midpoint correction were preferred.
% -------------------------------------------------------------------------

arguments
    obj 
    q1LOS
    medianLOS
    q3LOS
    K=60; %Default value
end

    [~, ~, time, A] = getPrediction(obj);
    A = A(:);                 % force admissions into a column vector
    T = length(A);

    % ---- 1) Convert empirical summaries to lognormal parameters ----
    %
    % If LOS ~ Lognormal(mu, sigma^2), then:
    %   mu    = log(median)
    %   sigma = (log(Q3) - log(Q1)) / (2 * 0.67449)
    %
    % where Q1 and Q3 are the 25th and 75th percentiles.
    [mu, sigma] = fitLognormalFromMedianIQR(medianLOS, q1LOS, q3LOS);

    % ---- 2) Build the discrete survival function p(k) ----
    %
    % p(k+1) corresponds to P(LOS >= k), for k = 0,...,K-1.
    %
    % Interpretation:
    %   p(1) = P(LOS >= 0) = 1
    %   p(2) = P(LOS >= 1)
    %   p(3) = P(LOS >= 2)
    %   ...
    %
    % This means p(k+1) is the fraction of patients expected to still be in
    % hospital k days after admission.
    p = makeDiscreteLognormal(mu, sigma, K);

    % ---- 3) Convolve admissions with survival probabilities ----
    %
    % prevalence(t) = sum over earlier admission cohorts that are still
    % present on day t.
    %
    % Using convolution gives the same result more compactly than a double
    % loop.
    prevalence = conv(A, p, 'full');
    prevalence = prevalence(1:T);

    % NOTE:
    % Because the lognormal distribution has no hard upper bound, K should
    % be chosen large enough so that the remaining tail probability beyond K
    % days is very small. In practice, it is useful to check p(end). If it
    % is still appreciable, K may need to be increased.
end


%% ------------------------------------------------------------------------
function [mu, sigma] = fitLognormalFromMedianIQR(medianLOS, q1LOS, q3LOS)
% fitLognormalFromMedianIQR
% -------------------------------------------------------------------------
% Fit lognormal parameters from empirical median and IQR.
%
% If X ~ Lognormal(mu, sigma^2), then:
%
%   mu = log(median)
%
% and
%
%   sigma = (log(Q3) - log(Q1)) / (2 * 0.67449)
%
% because the 25th and 75th percentiles of a standard normal distribution
% occur at approximately ±0.67449 standard deviations from the mean.
% -------------------------------------------------------------------------

    if medianLOS <= 0 || q1LOS <= 0 || q3LOS <= 0
        error('medianLOS, q1LOS, and q3LOS must all be positive.');
    end

    if ~(q1LOS <= medianLOS && medianLOS <= q3LOS)
        error('Inputs must satisfy q1LOS <= medianLOS <= q3LOS.');
    end

    mu = log(medianLOS);
    sigma = (log(q3LOS) - log(q1LOS)) / (2 * 0.67449);

    if sigma <= 0
        error('Computed sigma must be positive.');
    end
end


%% ------------------------------------------------------------------------
function p = makeDiscreteLognormal(mu, sigma, K)
% makeDiscreteLognormal
% -------------------------------------------------------------------------
% Construct a discrete survival function from a continuous lognormal LOS.
%
% OUTPUT:
%   p(k+1) = P(LOS >= k), for k = 0,...,K-1
%
% WHY THIS IS NEEDED:
%   The lognormal LOS model is continuous in time, but hospital admissions
%   and hospital prevalence are usually available as daily counts. Therefore,
%   we convert the continuous survival function into a daily survival vector.
%
% CURRENT CHOICE:
%   For k >= 1, we use
%
%       p(k+1) = 1 - logncdf(k, mu, sigma)
%
%   which corresponds to:
%
%       p(k) = P(LOS >= k)
%
%   on the integer-day grid.
%
% INTERPRETATION:
%   - k = 0  -> patient is certainly still present at the moment of admission
%   - k = 1  -> probability of still being hospitalized after 1 day
%   - k = 2  -> probability of still being hospitalized after 2 days
%   - etc.
%
% IMPORTANT SUBTLETY:
%   This is one reasonable discretization, but not the only one. Depending
%   on how daily prevalence is defined, one could instead use a midpoint
%   convention such as:
%
%       p(k+1) = 1 - logncdf(k + 0.5, mu, sigma)
%
%   or another day-binning convention.
%
%   So the code below is not "universally exact"; it is exact relative to
%   the chosen interpretation that p(k) = P(LOS >= k).
% -------------------------------------------------------------------------

    if sigma <= 0
        error('Lognormal distribution requires sigma > 0.');
    end

    p = zeros(K,1);

    for k = 0:(K-1)
        if k == 0
            % Everyone is still hospitalized at time 0, immediately after
            % admission.
            p(k+1) = 1;
        else
            % Survival probability k days after admission under the chosen
            % discretization convention.
            p(k+1) = 1 - logncdf(k, mu, sigma);
        end
    end
end