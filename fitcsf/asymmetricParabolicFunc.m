function CS = asymmetricParabolicFunc( p, SF )
% A 4-parameter asymmetric parabolic function
%
% Syntax:
%   CS = asymmetricParabolicFunc( p, SF )
%
% Description:
%   Implements the asymmetric parabolic function described in:
%
%       Chung ST, Legge GE. Comparing the shape of contrast sensitivity
%       functions for normal and low vision. Investigative ophthalmology &
%       visual science. 2016 Jan 1;57(1):198-207.
%
%   The function takes and returns values of spatial frequency and contrast
%   sensitivity in linear units. These values are log-transformed
%   internally for fitting. The parameters are:
%     - CSp: The peak contrast sensitivity value
%     - SFp: The spatial frequency at which CSp occurs
%     - widthL and widthR: The curvatures of the left and right branches
%           of the asymmetric parabolic function, respectively
%
% Inputs:
%   p                     - 1x4 vector. Parameters for the fit.
%   SF                    - 1xn vector. A vector of spatial frequencies
%                           in cycles / degree at which the function will
%                           be evaluated
%
% Outputs:
%   CS                    - 1xn vector. The contrast sensitivity values at
%                           the specified spatial frequencies.
%
% Examples:
%{
    p = [166, 2.5, 0.68, 1.28];
    SF = logspace(log10(0.2),log10(20),100);
	CS = asymmetricParabolicFunc(p,SF);
    loglog(SF,CS)
    xlabel('Spatial frequency [log cpd]');
    ylabel('Log contrast sensitivity');
%}

% Unpack the parameters and log-transform as needed
CSp = log10(p(1));
SFp = log10(p(2));
widthL = p(3);
widthR = p(4);

% Set the exponent to create parabolic functions
n = 2;

% Log-transform the spatial frequency
SF = log10(SF);

% Generate the contrast sensitivity values
CS = zeros(size(SF));
CS(SF<SFp) = CSp - (SF(SF<SFp)-SFp).^n * widthL^n;
CS(SF>=SFp) = CSp - (SF(SF>=SFp)-SFp).^n * widthR^n;

% Return the CS values to linear units and return
CS = 10.^CS;

end