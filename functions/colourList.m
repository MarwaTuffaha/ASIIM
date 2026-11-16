function [c] = colourList(n, colourSet)
% colourList Generates a list of distinguishable colors.
%
% Syntax:
%   c = colourList(n)
%   c = colourList(n, colourSet)
%
% Description:
%   The function `colourList` returns an array of colors for plotting purposes.
%   The number of colors returned is specified by the input `n`. The colors 
%   are chosen from a predefined set, which can be specified using the optional 
%   argument `colourSet`. If `colourSet` is 'large', the function returns a 
%   broader range of distinguishable colors.
%
% Input Arguments:
%   n         - Number of colors to return (scalar).
%   colourSet - Set of colors to use ('default' or 'large') (string, optional).
%               Default is 'large'.
%
% Output Arguments:
%   c - An n-by-3 matrix of RGB color values.
%
% Examples:
%   % Get 3 colors from the default set
%   colors = colourList(3);
%
%   % Get 5 colors from the 'large' set
%   colors = colourList(5, 'large');

arguments
    n
    colourSet = 'large';
end

if strcmp(colourSet, 'large')
    colour = distinguishable_colors(17,[0 0 0;1 1 1]);
else
    colour = [[0, 0.4470, 0.7410];
        [0.8500, 0.3250, 0.0980];
        [0.9290, 0.6940, 0.1250];
        [0.4940, 0.1840, 0.5560];
        [0.4660, 0.6740, 0.1880];
        [0.3010, 0.7450, 0.9330];
        [0.6350, 0.0780, 0.1840]];
end
c = colour(n,:);
end


