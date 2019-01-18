function [C,rc,modval,modls]=getdecodepth_Langlais(spec,ls,rplanet,showit)
% [C,rc,modval,modls]=getdecodepth_Langlais(spec,ls,rplanet,showit)
%
% SPEC MUST ONLY CONTAIN NON-ZONAL COEFFICIENTS
%
% Calculates the decorrelation depth described by Langlais et al (2014),
% eq. (5b) by doing a linear least squares fit to the logarithms.
% 
% INPUT:
%
% spec      spectrum
% ls        spherical-harmonic degrees for which the spectral values are
%           given. Must not contain 0!!! (log of zero is not a good thing)
% rplanet   radius of the planet
% showit    would you like to see a plot showing the fit? 1 for yes,
%           0 for no (default=0)
%
% OUTPUT:
%
% C         constant of flat spectrum
% rc        radial position for a flat spectrum 
%           (decorellation radial position)
% modlval   spectrum from fitting
% modls     degrees for fitted spectrum
% 
% Last modified by plattner-at-alumni.ethz.ch, 12/17/2018
  
if ls(1)==0
   error('Can not have degree 0')
end

if nargin<4
  showit=0;
end

spec=spec(:);
ls=ls(:);

x=2*ls+4;
y=-log(spec);

[m,q]=fitlin(x,y);



if showit
    plot(x,y,'x')
    hold on
    plot([x(1) x(end)],m*[x(1) x(end)]+q,'r--')
    hold off
end

C=exp(-q);
rc=rplanet/exp(m);

modls=ls;
modval=C*(rc/rplanet).^(2*ls+4);

end

function [m,q]=fitlin(x,y)
% [m,q]=fitlin(x,y)
%
% The simplest possible linear least squares

x=x(:);
y=y(:);

A=[x ones(size(x))];

v=(A'*A)\(A'*y);
m=v(1);
q=v(2);

end
