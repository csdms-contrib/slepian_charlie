function varargout=wattsandmoore(EL,tip,convo,wit,norma,r,ifn,mods,plans,degres)
% [dlola,degres,EL]=wattsandmoore(EL,tip,convo,wit,norma,r,ifn,mods,plans,degres) 
%
% Reproduces Figure 1 and Figure 3 from Watts and Moore, doi:10.1002/2017JB014571,
% i.e. a degree variance plot and map of bandlimited spatial expansion of
% the EGM2008 gravitational potential model by Pavlis et al., doi:10.1029/2011JB008916
%
% INPUT:
%
% EL      bandwidth (1 number, >=2, defaulted) or passband (2 numbers), defaulted
% tip     1 spatial map from PLOTPLM/PLM2XYZ
%         2 spectral rendition from PLM2SPEC [default]
% convo     1 free-air gravity in standard units
%         1e5 free-air gravity in mgal [default]
% wit     Option passed to PLM2POT
%           'nothing' when no reference is subtracted [default]
%           'WGS84' if your reference is WGS84 (even-degree zonals)
%           'spherical' only ever takes out the (0,0) coefficients
% norma   Option passed to PLM2SPEC for the sums of squares normalization
%          1 multiplication by (l+1) 
%          2 division by (2*l+1)
%          3 none, i.e. a scaling factor of 1 [default]
% r       Radius at which this is being evaluated, by default at the reference
% ifn     0 plots the SIGNAL strength under option tip=2 [default]
%         2 plots the NOISE strength under option tip=2
% mods    'EGM2008' as actually in the Watts and Moore paper [default]
%         'GMM3' as in a future paper
% plans   'Earth' as actually in the Watts and Moore paper [default]
%         'Mars' as in a future paper
% degres  Longitude/ latitude spacing, in degrees (for plotting maps)
%
% OUTPUT:
%
% dlola   Spatial expansion of the data, if you request this, you don't get
%         the figure, if tip===1 OR: 
%         Spectral density of the plot... if tip==2
% degres  Longitude/ latitude spacing, in degrees, OR:
%         the spherical harmonic degrees... if tip==2
% EL      Regurgitate the input
%
% SEE ALSO:
%
% PLM2POT, PLM2SPEC, PLOTPLM, PLM2XYZ
%
% EXAMPLE: 
%
%% Reproduces a piece of 10.1002/2017JB014571, Figure 1
% wattsandmoore([3 400],2,1e5,'nothing',3,[],0,'EGM2008','Earth',[])
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 03/18/2020

% Should figure out to keep degres constant or variable... for movie-type plots

% Default values
defval('EL',[3 400]);
defval('tip',2)
defval('convo',1e5)
defval('wit','nothing');
defval('norma',3);
defval('ifn',0);
% Model specification defaults
defval('mods','EGM2008');
defval('plans','Earth');

% Default pixel resolution - if you set to empty will get PLM2XYZ's default
defval('degres',[]);

% Unit awareness, not complete
switch convo
 case 1
  units='m/s^2';
 case 1e5
  units='mgal';
end

% Minimum spherical harmonic degree assumed in the file
minel=2;

% Those are the degrees requested
if size(EL)==1
  % All the possible degrees; the data starts at 2
  EL=[minel EL];
else
  % Adjust to take into account you must have minimum 2
  EL=[max(minel,EL(1)) EL(2)];
end

% Talk
disp(sprintf('\nWanting degrees %i to %i',EL(1),EL(2)))

% Create output filename for future use
fnpl=fullfile(getenv('IFILES'),'GRAVITY',mods,...
	      sprintf('%s_FreeAir_%i_%i.mat',mods,EL));

% Load, convert, expand to space, and save it or load it
if exist(fnpl,'file')~=2 
  % Get the potential coefficients from file
  if strcmp(mods,'EGM2008')
    % Additional model specification...
    egm=fralmanac(sprintf('%s_ZeroTide',mods),'SHM');
  else
    % No other funny business
    egm=fralmanac(sprintf('%s',mods),'SHM');
  end
    % Check the assumption that the data started only at minel
  diferm(minel,egm(1))

  % Get the reference GM product
  GM=fralmanac(sprintf('GM_%s',mods),plans);
  % Get the reference radius
  a=fralmanac(sprintf('a_%s',mods),plans);
  % Get the evaluation radius
  defval('r',a);
  
  % So we rescale all the coefficients by force to 2=free-air anomaly
  % Returns free-air gravity at r=a, as per readme.egm96, and Blakely Eq (7.2)
  egm=plm2pot(egm,r,GM,a,2,wit);
  % Unit conversion
  egm(:,3:end)=egm(:,3:end)*convo;

  % Remember that PLM2POT did extra work in supplying the 0 and 1 degrees
  % Restrict to the range you want, given where the data start
  egm=egm(addmup(EL(1)-1)+1-addmup(egm(1)-1):addmup(EL(end))-addmup(egm(1)-1),:);

  % Talk!
  disp(sprintf('\nGetting degree  %i order %i to degree %i order %i',...
	       egm(1,1:2),egm(end,1:2)))

  % Expand into spatial coordinates at appropriate degree resolution
  [dlola,lon,lat,~,degres]=plm2xyz(egm(:,[1:4]),degres);
  % Save all this information for subsequent retrieval
  save(fnpl,'dlola','egm','lon','lat','degres')
else
  disp(sprintf('Use preloaded file %s',fnpl))
  % Load what was has been calculated for fast access
  switch tip
   case 1
    % Only the expanded field
    load(fnpl,'dlola','degres')
    % If you don't make figure, dlola etc are returned
    % Create labels for future use
   case 2
    % Only the coefficients
    load(fnpl,'egm')
    % Talk!
    disp(sprintf('\nGetting degree  %i order %i to degree %i order %i',...
		 egm(1,1:2),egm(end,1:2)))
    % Spectral calculation of signal or noise - watch the normalization
    [sdl,l,bta,lfit,logy,logpm]=plm2spec(egm(...
        addmup(EL(1)-1)+1-addmup(egm(1)-1):addmup(EL(end))-addmup(egm(1)-1),...
        [1:4]+[0 0 ifn ifn]),norma);
    % If we don't make figure, return spectrum under the fake name
    dlola=sdl;
    degres=EL(1):EL(end);
  end
end

% Only make a figure if you don't request output
if ~nargout
  % Begin figure and create axis
  clf; ah=gca;
  
  % Now decide what figure to make
  switch tip
   case 1
    % Spatial map

    % Plot on the sphere 
    fig2print(gcf,'portrait')
    [r,c,ph]=plotplm(setnans(dlola),[],[],1,degres);
    
    % No continents on Mars
    if strcmp(plans,'Mars')
      delete([c{1} ph])
    else
      % But plot the hemispherical dichotomy
    end
    
    % Create labels for future use
    xxlabs=sprintf('%s free-air gravity anomaly in %s',mods,units);
    xlb=sprintf('spherical harmonic degrees %3.3i-%3.3i',EL);
    
    % Color scale
    kelicol

    % A relative scale that is a feast to the eyes
    caxis(round(10.^max(halverange(log(abs(r)),15)))*[-1 1]*1e5/convo)
    % An absolute scale that is a relative feast to the eyes
    caxis([-1 1]*1e2*1e5/convo)
    if strcmp(plans,'Mars')
      caxis([-1.5 1.5]*1e2*1e5/convo)
    end
    
        % Color bar
    cb=colorbar('hor');
    
    % Version control
    try
      axes(cb); 
      xt=title(xxlabs);
      xlabel(xlb)
      movev(xt,range(get(cb,'ylim'))/2)
      longticks(cb)
      movev(cb,-.15)
    catch
      cb.Label.String=sprintf('%s\n%s',xxlabs,xlb);
      cb.TickDirection='out';
      movev(cb,-.075)
    end
    
    shrink(cb,2,2)
    
    % Make it bigger to get a good bounding box
    set(ah,'camerav',6.5)
    movev([ah cb],.05)
    
    % Output to PDF
    figdisp([],sprintf('%i_%3.3i_%3.3i',tip,EL),[],2)

   case 2
    % Spectral plot

    fig2print(gcf,'portrait')

    % The power spectral density
    a=loglog(l,sdl,'o');
    hold on
    % The loglinear fit
    b=loglog(lfit,logy,'k-');
    hold off

    % Take that fit off here
    delete(b)
    
    % Create labels for future use
    xlabs='spherical harmonic degree';
    xxlabs='equivalent wavelength (km)';
    ylabs=sprintf('%s power spectral density [%s**2]',mods,units);

    % Cosmetics to match James' plot
    set(a,'MarkerFaceColor','k','MarkerSize',3,'MarkerEdgeColor','k')
    
    xlim(EL+[-0.5 40])
    longticks(gca)
    shrink(gca,1.333,1.075)

    % This needs to be data-dependent
    %   ylim([0.07 1e4]*1e5/convo)

    % The reference degrees you want plotted also
    nn=[12 33 400];
    % The reference degrees you wanted as well
    hold on
    for index=1:length(nn)
      pn(index)=plot([nn(index) nn(index)],ylim,'k--');
    end
    hold off

    % Labels
    ylabel(ylabs)
    xlabel(xlabs)

    % The degrees you want labeled
    if egm(1)<10 & egm(1,:)>100
      % Replace labels with meaningful ones
      els=[10 100]; ell={'10' '100'};
      set(ah,'XtickLabel',ell)
    end
    
    % Extra axis in equivalent wavelengths
    nlt=[10000 1000 100];
    [ax,xl,yl]=xtraxis(ah,round(jeans(nlt,0,1)),nlt,xxlabs);
    longticks(ax)
    % Output to PDF
    figdisp([],tip,[],2)
  end
end

% Create optional output
varns={dlola,degres,EL};
varargout=varns(1:nargout);
