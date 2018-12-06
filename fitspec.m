function [fspec,res]=fitspec(spec,ls,d,rplanet,rspec)
% [fspec,res]=fitspec(spec,ls,d,rplanet,rspec)
%
% Calculates a best fitting spectrum given a decorrelation depth
%
% INPUT:
%
% spec      Mauersberger-Lowes power spectrum for given degrees ls
% ls        degrees for the spectrum
% d         decorrelation radial position [typically close to rplanet]
% rplanet   planet radius
% rspec     radial position of spectrum
%
% OUTPUT:
%
% fspec     best fitting spectrum
% res       residual (spec - fspec)
%
% Last modified by plattner-at-alumni.ethz.ch, 4/9/2018

spec=spec(:);
ls=ls(:);


y=log(spec./(ls.*(ls+0.5).*(ls+1)));
x=ls;

m = 2*( log(d/rplanet) + log(rplanet/rspec) );
q = mean(y-m*x);

A=exp(q+2*log(d/rplanet)-4*log(rplanet/rspec));
fspec=A*ls.*(ls+0.5).*(ls+1).*(d/rplanet).^(2*ls-2).*(rplanet/rspec).^(2*ls+4);

res=spec-fspec;
