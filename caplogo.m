function [ah,ax]=caplogo(ahin,sord,loc,wid,hit,ang)
% [ah,ax]=caplogo(ahin,sord,loc,wid,hit,ang)
%
% Puts a thumbnail of a polar cap on another plot... some after-the-fact
% moving and shrinking may still be required
%
% INPUT:
%
% ahin      The axis handle on which you want the logo [default: gca] 
% sord      1 single cap [default]
%           2 double cap whereby the BACK of the bottom cap is shown
%           3 double cap whereby the FRONT of the bottom cap is shown
% loc       location, 'll', 'lr', 'ul', 'ur', 'um', 'ur' [default: 'lr']
%           OR two coordinates specifying the center
% wid       Width of the axis [default: 1/20]
% hit       Height of the axis [default: 1/20]
% ang       Angle of the cap [default: 20]
%
% OUTPUT:
%
% ah        The extra axis containing the top polar cap
% ax        The extra axis containing the bottom polar cap if available
%
% Last modified by fjsimons-at-alum.mit.edu, 2/19/2007

defval('ahin',gca)
defval('sord',1)
defval('loc','lr')
defval('wid',1/20)
defval('hit',1/20)

pp=getpos(ahin);
% See where the other axis is located
xl=[pp(1) pp(1)+pp(3)];
% Was wrong, had 1
yl=[pp(2) pp(2)+pp(4)];
% Axis width and height
rxl=pp(3);
ryl=pp(4);

xmrg=0; % Margin as ratio of xlim
ymrg=0; % Margin as ratio of ylim
%-----------------------------------------------

% Define mirroring operators
Mx=[ 1,-1];
My=[-1, 1];

% Define one box in upper right relative to midpoint
mid=boxmid([xl yl]);
% Define Left and Top corner, shifted by MID
X1=[xl(2)-rxl*xmrg-wid yl(2)-ryl*ymrg        ]-mid;
% Define Right and Bottom corner, shifted by MID
X2=[xl(2)-rxl*xmrg         yl(2)-ryl*ymrg-hit]-mid;

if isstr(loc)
  switch loc
   case 'ur'
    [B1,B2]=deal(X1,X2);
   case 'lr'
    [B1,B2]=deal(Mx.*X1,Mx.*X2);
   case 'ul'
    [B1,B2]=deal(My.*X1,My.*X2);
   case 'll'
    [B1,B2]=deal(My.*Mx.*X1,My.*Mx.*X2);
   case 'um'
    [B1,B2]=deal([-wid/2 X1(2)],[wid/2  X2(2)]);
   case 'lm'
    [B1,B2]=deal([-wid/2 -X1(2)],[wid/2  -X2(2)]);
  end
else
  B1=[loc(1)-wid/2 loc(2)+hit/2];
  B2=[loc(1)+wid/2 loc(2)-hit/2];
  mid=0;
end

% Shift the MID back in
B1=B1+mid;
B2=B2+mid;
bcor=[B1(1) B2(1) B1(2) B2(2)];
tcor=boxmid(bcor);

% Create axis % Used to be wrong, not sorted
ah=axes('position',lrtb2ext([sort(bcor(1:2)) bcor([4 3])]));

% Plot the double or the single polar cap
if sord==1
  % Plot in sphere etc
  [h,cord]=circ(1); delete(h)
  % Remove this if you want a transparent background
  ol2=fill3(-1*ones(size(cord(:,1))),cord(:,2),cord(:,1)-0.32,'w');
  set(ol2,'edgec','w')
  hold on
  ol=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
  [h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
  eq=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),'-');
  [h,cord]=circ(1,[pi/2 3*pi/2]); delete(h)
  eqd=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
  % Plot polar cap
  defval('ang',20);
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  [h,cord2]=circ(1,[ang 180-ang]*pi/180); delete(h)
  eqx=fill3([cord(:,1) ;  ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ;  ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang*pi/180) ; ... 
	     ; cord2(:,2)]',...
	    grey);
  hold off
  viewpars
  ax=NaN;
elseif sord==2
  [h,cord]=circ(1); delete(h)
  % Remove this if you want a transparent background
  ol2=fill3(-1*ones(size(cord(:,1))),cord(:,2),cord(:,1)-0.32,'w');
  set(ol2,'edgec','w')
  hold on
  ol=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
  [h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
  eq=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
  [h,cord]=circ(1,[0 3*pi/2]); delete(h)
  eqd=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
  % Plot polar cap on TOP
  defval('ang',20);
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  eq(2)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  % Make this grey
  [h,cord2]=circ(1,[ang 180-ang]*pi/180); delete(h)
  eqx=fill3([cord(:,1) ;  ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ;  ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang*pi/180) ; ... 
	     ; cord2(:,2)]',...
	    grey);
  viewpars
  % Plot BACK of the polar cap on BOTTOM on extra axis
  [ax,axl,loc]=laxis(ah,0);
  axes(ax)
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  hold on
  eq(3)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  [h,cord2]=circ(1,[ang 180-ang]*pi/180); delete(h)
  eqx=fill3([cord(:,1) ;  ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ;  ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang*pi/180) ; ... 
	     ; cord2(:,2)]',...
	    grey);
  set(eqx,'LineS','-')
  hold off 
  axis([get(ah,'xlim') get(ah,'ylim') get(ah,'zlim')])
  view(90,17.5)
  zlim([0.7 2.7])
  axis off
elseif sord==3
  [h,cord]=circ(1); delete(h)
  % Remove this if you want a transparent background
  ol2=fill3(-1*ones(size(cord(:,1))),cord(:,2),cord(:,1)-0.32,'w');
  set(ol2,'edgec','w')
  hold on
  ol=plot3(zeros(size(cord(:,1))),cord(:,2),cord(:,1));
  [h,cord]=circ(1,[-pi/2 pi/2]); delete(h)
  eq=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))));
  [h,cord]=circ(1,[0 3*pi/2]); delete(h)
  eqd=plot3(cord(:,1),cord(:,2),zeros(size(cord(:,1))),':');
  % Plot polar cap on TOP
  defval('ang',20);
  % Make an equatorial half circle
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  % But plot it up
  eq(2)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  % Define another boundary
  [h,cord2]=circ(1,[ang 180-ang]*pi/180); delete(h)
  % And fill the enclosed space with grey
  eqx=fill3([cord(:,1) ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang*pi/180) ; cord2(:,2)]',... 
	    grey);
  viewpars
  % Plot bottom polar cap on TOP
  ang=-ang;
  % Make an equatorial half circle
  [h,cord]=circ(cos(ang*pi/180),[-pi/2 pi/2]); delete(h)
  % But plot it at the right height
  eq(3)=plot3(cord(:,1),cord(:,2),ones(size(cord(:,1)))*sin(ang*pi/180));
  % Define another boundary
  [h,cord2]=circ(1,[-ang 180+ang]*pi/180); delete(h)
  % Make this grey
  eqx=fill3([cord(:,1) ; zeros(size(cord2(:,1)))]',...
	    [cord(:,2) ; cord2(:,1)]',...
	    [ones(size(cord(:,1)))*sin(ang*pi/180) ; -cord2(:,2)]',...
	    grey);
  hold off 
  ax=NaN;
else
  error('Specify valid option')
end

set(ol,'color','k')
set(eq,'color','k','lines',':')
set(eqd,'color','k','lines',':')

% Maybe? Delete the white disk?
%delete(ol2)
% Delete the equators
delete([eq eqd])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewpars
% Set viewing parameters
% Really would need another rotation around x
%view(90,17.5); axis equal
view(90,17.5); axis equal
axis([-1.01 1.01 -1.01 1.01 -1.01 1.01]); axis off
