function slola=wattsandmoore(EL,fig,convo,wit,norma,r,ifn)
% slola=wattsandmoore(EL,fig,convo,wit,norma,r,ifn) 
%
% Reproduces Figure 1 and Figure 3 from Watts and Moore, doi:10.1002/2017JB014571,
% i.e. a degree variance plot and map of bandlimited spatial expansion of
% the EGM2008 gravitational potential model by Pavlis et al., doi:10.1029/2011JB008916
%
% INPUT:
%
% EL      bandwidth (1 number, >=2, defaulted) or passband (2 numbers), defaulted
% fig     1 spatial map from PLOTPLM/PLM2XYZ
%         2 spectral rendition from PLM2SPEC [default]
% convo   1 free-air gravity in standard units
%         2 free-air gravity in mgal [default]
% wit     Option passed to PLM2POT
%           'nothing' when no reference is subtracted [default]
%           'WGS84' if your reference is WGS84 (even-degree zonals)
%           'spherical' only ever takes out the (0,0) coefficients
% norma   Option passed to PLM2SPEC for the sums of squares normalization
%          1 multiplication by (l+1) 
%          2 division by (2*l+1)
%          3 none, i.e. a scaling factor of 1 [default]
% r       Radius at which this is being evaluated, by default at the reference
% ifn     0 plots the SIGNAL strength under option fig=2 [default]
%         2 plots the NOISE strength under option fig=2
%
% OUTPUT:
%
% slola   Spatial expansion of the data
%
% SEE ALSO PLM2POT, PLM2SPEC, PLOTPLM, PLM2XYZ
%
% Tested on 8.3.0.532 (R2014a) and 9.0.0.341360 (R2016a)
% Last modified by fjsimons-at-alum.mit.edu, 12/6/2018

% Should figure out to keep degres constant or variable... for movie-type plots

% Default values
defval('EL',400);
defval('fig',2)
defval('convo',1e5)
defval('wit','nothing');
defval('norma',3);
defval('ifn',0);

% Default pixel resolution
defval('degres',0.25);

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

% Create filename for future use
fnpl=fullfile(getenv('IFILES'),'GRAVITY',...
	      sprintf('EGM2008_FreeAir_%i_%i.mat',EL));

% Load, convert, expand to space, and save it or load it
if exist(fnpl,'file')~=2
  % Get the potential coefficients from file
  egm=fralmanac('EGM2008_ZeroTide','SHM');
  % Check the assumption that the data started only at minel
  diferm(minel,egm(1))

  % Get the reference GM product
  GM=fralmanac('GM_EGM2008','Earth');
  % Get the reference radius
  a=fralmanac('a_EGM2008','Earth');
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
  [slola,lon,lat]=plm2xyz(egm(:,[1:4]),degres);
  % Save all this information for subsequent retrieval
  save(fnpl,'slola','egm','lon','lat')
else
  disp(sprintf('Use preloaded file %s',fnpl))
  % Load what was has been calculated for fast access
  switch fig
   case 1
    % Only the expanded field
    load(fnpl,'slola')
    % Create labels for future use

   case 2
    % Only the coefficients
    load(fnpl,'egm')
    % Talk!
    disp(sprintf('\nGetting degree  %i order %i to degree %i order %i',...
		 egm(1,1:2),egm(end,1:2)))
  end
end

 % Begin figure and create axis
 figure(1); clf; ah=gca;

 % Now decide what figure to make
 switch fig
  case 1
   % Spatial map

   % Plot on the sphere 
   fig2print(gcf,'portrait')
   [r,c,ph]=plotplm(setnans(slola),[],[],1,degres);

   % Create labels for future use
   xxlabs=sprintf('EGM2008 free-air gravity anomaly in %s',units);
   xlb=sprintf('spherical harmonic degrees %3.3i-%3.3i',EL);
   
   % Color scale
   kelicol

   % A relative scale that is a feast to the eyes
   caxis(round(10.^max(halverange(log(abs(r)),15)))*[-1 1]*1e5/convo)
   % An absolute scale that is a feast to the eyes
   caxis([-1 1]*1e2*1e5/convo)

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
        
   figdisp([],sprintf('%i_%3.3i_%3.3i',fig,EL),[],2)

  case 2
   % Spectral plot

   % Spectral calculation of signal or noise - watch the normalization
   [sdl,l,bta,lfit,logy,logpm]=plm2spec(egm(:,[1:4]+[0 0 ifn ifn]),norma);

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
   ylabs=sprintf('EGM2008 power spectral density in [%s]^2',units);

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
   figdisp([],fig,[],2)
 end
 
