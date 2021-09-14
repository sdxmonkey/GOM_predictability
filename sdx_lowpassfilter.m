% LOWPASSFILTER - Constructs a low-pass butterworth filter.
%
% usage: f = sdx_lowpassfilter(sze, dx, cutoff, n)
% 
% where: sze    is a two element vector specifying the size of filter 
%               to construct.
%        dx     is the spatial resoulution of data
%        cutoff is the cutoff frequency of the filter 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%
% The frequency origin of the returned filter is at the corners.
%
% See also: HIGHPASSFILTER, HIGHBOOSTFILTER, BANDPASSFILTER
%

% Peter Kovesi   pk@cs.uwa.edu.au
% Department of Computer Science & Software Engineering
% The University of Western Australia
%
% October 1999
%
% Modified Rob Gaddi   gaddi@rice.edu
% ELEC 301
% Rice University
%
% December 2001
%
% Modified by sdxmonkey sdxmonkey@gmail.com to produce filter on
% non-normorlozed frequecy space
%

function f = sdx_lowpassfilter(sze, dx, cutoff, n)
    
    if cutoff < 0
	error('cutoff frequency must be positive');
    end
    
    if rem(n,1) ~= 0 || n < 1
	error('n must be an integer >= 1');
    end
    
    %  Modification ELEC 301 Project Group, Dec 2001
    %  Original code [rows, cols] = sze was not accepted by Matlab
    rows = sze(1);
    cols = sze(2);
    %  End Alteration

    % X and Y matrices with ranges normalised to +/- 0.5
    x =  (ones(rows,1) * [1:cols]  - (fix(cols/2)+1))/cols;
    x = x/dx;
    y =  ([1:rows]' * ones(1,cols) - (fix(rows/2)+1))/rows;
    y = y/dx;
    
    radius = sqrt(x.^2 + y.^2);        % A matrix with every pixel = radius relative to centre.
    
    %  Alteration, ELEC 301 Project Group, Dec 2001
    %  Original code fftshifted the filter before output.  Since
    %  imFFT and imIFFT already shift, the output should remain low-centered.
    f = 1 ./ (1.0 + (radius ./ cutoff).^(2*n));   % The filter
    %  End Alteration