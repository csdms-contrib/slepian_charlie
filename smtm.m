function [sdlMT,l]=smtm(L,R,dlola,degres)
% [sdlMT,l]=smtm(L,R,dlola,degres)
%
% Bandlimited, geographically localized, spherical multitaper
% power-spectral density estimate of a spatially gridded global data set,
% as in Dahlen and Simons, doi: 10.1111/j.1365-246X.2008.03854.x
%
% INPUT:
%
% L         multitaper bandwidth (keep it small)
% R         region of interest, e.g., 'alloceans', 'contshelves', 'africa', ...
% dlola     global spatial expansion, on a complete longitude and latitude
%           grid, as coming out of PLM2XYZ or some such function (see example).
% degres    the grid step in longitudinal degrees
%
% OUTPUT:
%
% sdlMT      multitaper spectral density estimate
%
%% Last modified by fjsimons-at-alum.mit.edu, 12/7/2018

% Default values
defval('R','contshelves')
defval('L',8)
defval('J',(L+1)^2)
defval('xver',1)
defval('degres',1)

% The example data are the EGM2008 zero-tide non-WGS84 corrected free-air
% gravity field at the reference radius for the EGM2008 model
defval('dlola',[])
defval('degres',[])
if isempty(dlola) & isempty(degres)
  [dlola,degres]=wattsandmoore([2 400]);
end

% This is the sum of the eigenvalues
N=spharea(R)*(L+1)^2;

% Get the spherical harmonic expansion coefficients of the Slepian basis
% functions in one big matrix
%[Glma,V,EL,EM,N]=glmalpha(R,L,[],[],[],[],J,[],1);
% Might be issues, try LOCALIZATION
[V,C,dels,dems,XY,Klmlmp,Glma]=localization(L,R,[],J,1,0);

% Take a look using PLOTSLEP, CONTSHELVES, see also LOCALIZATION
if xver==1
  for index=1:(L+1)^2
    clf
    % This should return the spatial taper
    % gialpha=plotslep(Glma,index);
    gialpha{index}=plotplm([dels dems C{index}],[],[],4,degres); hold on
    % This plots the region, for now, hardcoded
    contshelves([],1)

    % Labels
    title(sprintf('%s = %i ; N = %i out of %i','\alpha',...
		  index,round(N),(L+1)^2))
    pause
  end
  keyboard
end


% Prepare for arrival
%sdl=nan(egm(end,1),length(V));

% For all them
for index=1:length(V)
    disp(index)
    % Taper the data also fake it
    tapdat=dlola.*kindeks(gialpha{index},1:size(dlola,2));
    % imagesc(tapdat);
    % Spherical transform
    lmcosi(:,:,index)=xyz2plm(tapdat,max(EL));
    % Normalized sum of squares
    norma=3
    [sdl(:,index),l]=plm2spec(lmcosi(:,:,index),norma);
    loglog(l,sdl(:,index),'+');         
    drawnow
end

% At the end, average the result
sdlMT=mean(sdl,2)/sum(V);
