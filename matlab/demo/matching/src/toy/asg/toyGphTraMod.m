function [mess, Varss] = toyGphTraMod(shp, k, nF)
% Generate latent graph trajectory.
%
% Input
%   shp     -  graph shape, 1 | ...
%              1: random graph
%   k       -  #nodes
%   nF      -  #frames
%
% Output
%   mess    -  mean of Gaussian, 1 x nF (cell), 1 x k (cell), 2 x 1
%   Varss   -  variance of Gaussian, 1 x nF (cell), 1 x k (cell), 2 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% random graph
if shp == 1
    % parameter
    wei = .7;
    
    % Gauss for each node
    [mess, Varss] = cellss(1, nF);
    
    % first frame
    [mess{1}, Varss{1}] = modGausses(k, ...
                                     'dis', 'unif', ...
                                     'maMiD', 'n', ...
                                     'rBd', [.007, .02]);
    % per frame
    for iF = 2 : nF
        % mean
        mess{iF} = cell(1, k);
        
        % per node
        for c = 1 : k
            % velocity from last frame
            if iF == 2
                dx0 = [0; 0];
            else
                dx0 = mess{iF - 1}{c} - mess{iF - 2}{c};
            end
            
            % current velocity
            dx = randn(2, 1);

            % combine
            mess{iF}{c} = mess{iF - 1}{c} + wei * dx0 + (1 - wei) * dx;
        end

        % varaiance
        [~, Varss{iF}] = modGausses(k, ...
                                    'dis', 'unif', ...
                                    'maMiD', 'n', ...
                                    'rBd', [.007, .02]);
    end
    
else
    error('unknown shape: %d', shp);
end
