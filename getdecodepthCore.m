function [K,r,modval,modls]=getdecodepthCore(spec,ls,rplanet,showit)
% [K,r,modval,modls]=getdecodepthCore(spec,ls,rplanet,showit)
%
% THIS ONE IS USEFUL FOR THE CORE. THE OTHER ONE IS FOR CRUSTAL FIELDS
%  
% Calculates the decorrelation depth described by Voorhies et al (2002),
% eq. (16) by doing a linear least squares fit to the logarithms.
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
% K         power factor (see Voorhies et al paper)
% r         radial position of the random dipoles 
%           (decorellation radial position)
% modlval   spectrum from fitting
% modls     degrees for fitted spectrum
% 
% Last modified by plattner-at-alumni.ethz.ch, 12/18/2018
  
if ls(1)==0
   error('Can not have degree 0')
end

if nargin<4
    showit=0;
end

spec=spec(:);
ls=ls(:);

y=log(spec.*(ls+1).*ls./(ls+0.5));
x=2*ls+4;

[m,q]=fitlin(x,y);

if showit
    plot(x,y,'x')
    hold on
    plot([x(1) x(end)],m*[x(1) x(end)]+q,'r--')
end

r=exp(m)*rplanet;
K=exp(q);

modval=K*(ls+0.5)./(ls.*(ls+1)).*((r/rplanet).^(2*ls+4));
modls=ls;

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
