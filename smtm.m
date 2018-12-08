function [sdlMT,l,N]=smtm(L,R,dlola,degres,EL,xver)
% [sdlMT,l,N]=smtm(L,R,dlola,degres,xver)
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
% EL        data bandwidth (1 number, >=2, defaulted) or passband (2 numbers)
% norma     Option passed to PLM2SPEC for the sums of squares normalization
%           1 multiplication by (l+1) 
%           2 division by (2*l+1)
%           3 none, i.e. a scaling factor of 1 [default]
%
% OUTPUT:
%
% sdlMT      multitaper spectral density estimate
% l          the spherical harmonic degrees at which this is evaluated
% N          the rounded Shannon number
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)

%% Last modified by fjsimons-at-alum.mit.edu, 12/7/2018

% Default values
defval('R','contshelves')
defval('L',8)
defval('J',(L+1)^2)
defval('xver',1)
defval('norma',3)

% The example data are the EGM2008 zero-tide non-WGS84 corrected free-air
% gravity field at the reference radius for the EGM2008 model
defval('dlola',[])
defval('degres',[])
defval('EL',[])
if isempty(dlola) || isempty(degres) || isempty(EL)
  [dlola,degres,EL]=wattsandmoore([2 400],1);
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
    % bob=0;
    clf
    % This should return the spatial taper
    % gialpha=plotslep(Glma,index);
    gialpha{index}=plotplm([dels dems C{index}],[],[],4,degres); hold on
    % bob=bob+gialpha{index}.^2*V(index);
    % This plots the region, for now, hardcoded
    contshelves([],1)

    % Labels
    title(sprintf('%s = %i ; %s = %5.3f ; N = %i out of %i','\alpha',...
		  index,'\lambda',V(index),round(N),(L+1)^2))
    %pause(0.6)
  end
end

% imagefnan([0 90],[360 -90],bob)

% Prepare for arrival
howmany=round(N);
sdl=nan(max(EL)+1,howmany);
l=0:max(EL);

clf

% Round the Shannon number for output
N=round(N);

% For all them
parfor index=1:size(sdl,2)
    % disp(index)
    % Taper the data 
    tapdat=dlola.*gialpha{index};
    % imagesc(tapdat);
    % Spherical transform
    lmcosi(:,:,index)=xyz2plm(tapdat,max(EL));
    % Normalized sum of squares time the eigenvalues for averaging later
    sdl(:,index)=V(index)*plm2spec(lmcosi(:,:,index),norma);
    % Take a quick look without the eigenvalue
    % loglog(l,sdl(:,index)/V(index),'+');         
    % hold on
    % drawnow
end

% At the end, AVERAGE the result
sdlMT=size(sdl,2)*nanmean(sdl,2)/sum(V(1:size(sdl,2)));

% And make a final plot
loglog(l,sdlMT,'ko','MarkerF','k');

% And then we need to produce the variance estimates
% mtvar
keyboard
