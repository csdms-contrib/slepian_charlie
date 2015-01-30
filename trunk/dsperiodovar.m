function dsperiodovar(ELS)
% DSPERIODOVAR(ELS)
% Obsolete, but don't trash just yet

% Plots the ratio of cut-sky periodogram variance to whole-sky estimation
% variance for a variety of observable patches and in fuction of the
% degree, with asymptotics superimposed.
%
% INPUT:
%
% ELS      A (range of) spherical harmonic degree(s)
%
% Last modified by fjsimons-at-alum.mit.edu, 10/02/2006

% Size of the cut patch
defval('ELS',5:5:50);
THint=0.1;
TH=[1:THint:20 20:1:180];

% Find wavelength of best-fitting half box according to Jeans
thel=180./sqrt(ELS.*(ELS+1));

% Initialize arrays
v=repmat(NaN,length(ELS),length(TH));

% First time return the spectra as well
[v(1,:),v2,v3,v4,TH,A,Bl,bels]=periodovar(ELS(1),TH);
for index=2:length(ELS)
  v(index,:)=periodovar(ELS(index),TH,[],[],Bl,bels);
end

% % Find appropriate Lmax? 
% for in=1:size(Bl,2)
%   sk=Bl(:,in)<max(Bl(:,in))/1000;
%   if sum(sk)
%     LmaxN(in)=indeks(find(sk),1);
%   else
%     LmaxN(in)=NaN;
%   end
% end
% if any(Lmax<LmaxN)
%   error('Reconsider and make Lmax bigger')
% end

% Do the plotting
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(1)=subplot(211);
p1=plot(TH,v,'Color',grey);

hold on
p2=plot(TH,v2,'k-');
p3=plot(TH,v3,'k--');
xlim([1 15])
ylim([1 60])
set(ah(1),'ytick',[1 10:10:60])
xl(1)=xlabel(sprintf('colatitudinal radius %s','\Theta'));
yl(1)=ylabel('variance ratio');
set(ah(1),'xtick',[1 3:3:15])
aa=plot([thel' thel'],[1 60],'k:');
%lk=sort([1 thel 15]);
%set(ah(1),'xtick',lk)
%set(ah(1),'xtickl',round(lk*10)/10,'xgrid','on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(2)=subplot(212);
p4=plot(TH,v,'Color',grey);

hold on
p5=plot(TH,v2,'k-');
p6=plot(TH,v3,'k--');
p7=plot(TH,v4,'k--');
plot([0 180],[1 1],'k:')
plot([90 90],[0 5],'k:')
plot(90,1,'o','markerf','k','markere','k','markers',3)
plot(180,1,'o','markerf','k','markere','k','markers',3)
bb=plot([thel' thel'],[0.5 3.5],'k:');
hold off
xl(2)=xlabel(sprintf('cap radius %s','\Theta'));
yl(2)=ylabel('variance ratio');
set(ah(2),'xtick',[15:15:180])

% Cosmetic adjustments
xlim([15 180])
ylim([0.5 3.25])

fig2print(gcf,'portrait')
longticks(ah,2)
delete(xl(1))

set([p1; p4],'LineW',1)
serre(ah(1:2),1/2,'down')
deggies(ah(1:2),1)

figdisp



