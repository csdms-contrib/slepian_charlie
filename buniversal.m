function [Mro,lp]=buniversal(TH,Lmax,l,sord)
% [Mro,lp]=BUNIVERSAL(TH,Lmax,l,sord)
%
% Approximations to the boxcar coupling function, which are independent of
% the area of the concentration region.
%
% Dahlen & Simons (2008), Eq. 65, doi: 10.1111/j.1365-246X.2008.03854.x
%
% INPUT:
%
% TH     Colatitudinal radius of the cap, in degrees
%        Colatitudinal halfwidth of the belt, degress
% Lmax   Maximum degree of this non-bandlimited boxcar function
% l      The supposed degree index of the coupling matrix
% sord   1 Single cap of diameter 2TH [default]
%        2 Double cap left by subtracting belt of width 2TH
%        3 Equatorial belt of width 2TH
%
% OUTPUT:
%
% Mro    Coupling matrix values centered around the diagonal [crude approximation]
% lp     The supposed column index (for plotting purposes only)
%
% See also UNIVERSAL, for the same result of the multitaper estimate
%
% Last modified by fjsimons-at-alum.mit.edu, 03/26/2020

defval('TH',10)
defval('Lmax',100)
defval('l',0)
defval('sord',1)

Lpot=(Lmax+1)^2;

% Calculates the areas of these regions
A=4*pi*spharea(TH,sord);

% Get block-ordered Xlm functions
[X2,jk1,jk2,ems]=xlm(0:Lmax,[],pi/2,0,[],1);

% Get the power spectrum of the boxcar
[Bl,dels]=bpboxcap(TH,Lmax,[],0,sord);

% Now need to repeat this Lmax, Lmax-1, Lmax-2, and so on, ... times
% There may be a better way to do the next two lines
B2=gamini(Bl,2*(0:Lmax)+1);
[EM,EL,mz,blkm]=addmout(Lmax);
B2=B2(blkm)';

% Square, sum over all degrees and normalize
emo=indeks((repmat(0:Lmax,2,1)'*diag([-1 1]))','2:end');
Mro=sparse(gamini(1:(2*Lmax+1),Lmax-abs(emo)+1'),...
		     1:Lpot,1,2*Lmax+1,Lpot)*(B2.*X2.^2);
% Resort Mro to bring out the symmetry around the diagonal
[s,i]=sort(emo);
Mro=Mro(i)/A*4*pi;
lp=l-Lmax-1+[1:length(Mro)];



