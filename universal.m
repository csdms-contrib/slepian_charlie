function [Mro,lp,Mra,lpp]=universal(L,l)
% [Mro,lp,Mra,lpp]=UNIVERSAL(L,l)
%
% Approximations to the multitaper (with uniform weights) coupling function,
% independent of the area of the concentration region or the exact power
% spectrum of the Slepian tapers used in the spectral density estimate.
%
% Dahlen & Simons (2008), Eqs 145-146, doi: 10.1111/j.1365-246X.2008.03854.x
%
% INPUT:
%
% L      Bandwidth of the spectral windows used
% l      The supposed row index of the coupling matrix
%
% OUTPUT:
%
% Mro    Coupling matrix values centered around the diagonal [crude approximation]
% lp     The supposed column index (for plotting purposes only)
% Mra    Coupling matrix values centered around the diagonal [better approximation]
% lpp    The supposed column index (contained in the better approximation)
%
% Last modified by fjsimons-at-alum.mit.edu, 03/25/2020

% There's a slight redundancy here with the symmetry, but let's hold off
% optimizing this the moment. The crude approximation involves an
% asymptotic approximation of the better one, which merely replaces the
% sum of the power spectrum by the area. For the crude approximation, the
% only thing that matters is the distance to the diagonal and the row
% index serves merely for plotting purposes; the better approximation
% does depend on the row index l.

defval('L',18)
defval('l',0)

Lpot=(L+1)^2;

% Get block-ordered Xlm functions
[X2,jk1,jk2,ems]=xlm(0:L,[],pi/2,[],[],1);

% Square, sum over all degrees and normalize
emo=indeks((repmat(0:L,2,1)'*diag([-1 1]))','2:end');
Mro=4*pi/Lpot*sparse(gamini(1:(2*L+1),L-abs(emo)+1'),...
		     1:Lpot,1,2*L+1,Lpot)*X2.^2;
% Resort Mro to bring out the symmetry around the diagonal
[s,i]=sort(emo);
Mro=Mro(i);
lp=l-L-1+[1:length(Mro)];

% Now on to the approximation to which the previous is an approximation
if nargout>=3
  % This is the degree for the calculation
  lpp=max(0,l-L):l+L+1;
  Mra=repmat(NaN,length(lpp),1);
  for index=1:length(lpp)
    Mra(index)=(2*lpp(index)+1)*sum((2*[0:L]+1).*wigner0j(L,l,lpp(index)).^2);
  end
  % Check last one is indeed zero
  difer(Mra(end))
  % This one isn't symmetric around the diagonal
  Mra=Mra/Lpot;
  % But this needs to be corrected for the degree in the plotting scheme
  lpp=lpp-l;
end





