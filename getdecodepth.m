function [A,d,modval,modls]=getdecodepth(spec,ls,rplanet,rspec,showit)
% [A,d,modval,modls]=getdecodepth(spec,ls,rplanet,rspec,showit)
%
% Calculates the decorrelation depth described by Voorhies et al (2002),
% eq. (5b) by doing a linear least squares fit to the logarithms.
% 
% INPUT:
%
% spec      spectrum
% ls        spherical-harmonic degrees for which the spectral values are
%           given. Must not contain 0!!! (log of zero is not a good thing)
% rplanet   radius of the planet
% rspec     radial position at which the spectrum is given
% showit    would you like to see a plot showing the fit? 1 for yes,
%           0 for no (default=0)
%
% OUTPUT:
%
% A         power factor (see Voorhies et al paper)
% d         radial position of the random dipoles 
%           (decorellation radial position)
% modlval   spectrum from fitting
% modls     degrees for fitted spectrum
% 
% Last modified by plattner-at-alumni.ethz.ch, 01/17/2017
  
if ls(1)==0
   error('Can not have degree 0')
end

if length(nargin<5)
    showit=0;
end

spec=spec(:);
ls=ls(:);

y=log(spec./(ls.*(ls+0.5).*(ls+1)));
x=ls;

[m,q]=fitlin(x,y);

if showit
    plot(x,y,'x')
    hold on
    plot([x(1) x(end)],m*[x(1) x(end)]+q,'r--')
end

d=exp(0.5*(m-2*log(rplanet/rspec)))*rplanet;
A=exp(q+2*log(d/rplanet)-4*log(rplanet/rspec));

modval=A*ls.*(ls+0.5).*(ls+1).*(d/rplanet).^(2*ls-2).*(rplanet/rspec).^(2*ls+4);
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