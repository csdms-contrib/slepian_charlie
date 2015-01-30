function dsperiodovar2(ELS,Lmax)
% DSPERIODOVAR2(ELS,Lmax)
% Obsolete, but don't trash just yet

% Plots the ratio of cut-sky periodogram variance to whole-sky estimation
% variance for a variety of observable patches and in function of the
% degree, with asymptotics superimposed.
%
% INPUT:
%
% ELS      A (range of) spherical harmonic degree(s)
%
% Last modified by fjsimons-at-alum.mit.edu, 04/29/2007

% Size of the cut patch
defval('ELS',[0:2:50]);
defval('Lmax',100);
THint=2.5;
TH=[0.1 THint:THint:90];

% Find wavelength of best-fitting half box according to Jeans
warning off
thel=180./sqrt(ELS.*(ELS+1));
warning on
if ELS(1)==0
  thel(1)=90;
end

% Initialize arrays
v=repmat(NaN,length(ELS),length(TH));

% First time return the spectra as well
[v(1,:),v2,v3,v4,TH,A,Bl,bels]=periodovar(ELS(1),TH);
for index=2:length(ELS)
  v(index,:)=periodovar(ELS(index),TH,[],[],Bl,bels);
  v(index,:)=v(index,:);
end

% Find appropriate Lmax? 
for in=1:size(Bl,2)
  sk=Bl(:,in)<max(Bl(:,in))/1000;
  if sum(sk)
    LmaxN(in)=indeks(find(sk),1);
  else
    LmaxN(in)=NaN;
  end
end
if any(Lmax<LmaxN)
  error('Reconsider and make Lmax bigger')
end

% Do the plotting
clf

h=mesh(TH,ELS,v);
set(h,'edgec','k')

% Cosmetic adjustments
axis tight; xlim([0 90]); ylim([0 50]); zlim([1 100])
set(gca,'xgrid','off','ygrid','off','zgrid','off')
set(gca,'xtick',[TH(1) 10:10:90])
set(gca,'ztick',[1 10:10:100])
xl=xlabel(sprintf('colatitudinal radius %s','\Theta'));
yl=ylabel('spherical harmonic degree l');
zl=zlabel('variance ratio');
a=40;
e=15;
view(40,25)

set(xl,'rotat',-13)
set(yl,'rotat',21)

movev(xl,7.5); moveh(xl,-20)
moveh(yl,-10); 
box on

fig2print(gcf,'portrait')
longticks(gca)
deggies(gca,1)

figdisp

% What is the variance ratio at the box size appropriate for l?
% for index=1:length(ELS)
%   [vapp(index),vass(index)]=periodovar(ELS(index),thel(index));
% end
% What is vapp-vass? For various fractions of thel? Not very informative.

keyboard
