function varargout=localspectrum2(data,lon,lat,Lmax,Ltap,dom,Jmax,rotcoord,Nspec)
% [spec,specvar,spectap]=localspectrum2(data,lon,lat,Ltap,dom,Jmax,rotcoord,Nspec)
%
% Calculates the local multitaper spectrum using the formula described by
% Dahlen & Simons (2008), eq. (130) by multiplying the evaluated function
% with multitapers evaluated at the same locations as the function. 
%
% INPUT:
%
% data      Values of the function for which you want to calculate the
%           local spectrum
% lon       longitudinal locations of the function value points 
%           (0<=lon<=360)
% lat       latitudinal locations of the function value points 
%           (-90<=lat<=90)
% Lmax      Maximum spherical-harmonic degree for which you want to
%           calculate the regional spectral value
% Ltap      Maximum spherical-harmonic degree for the multitapers with
%           which you want to calculate the local spectrum. This defines
%           the trade-off between spectral and spatial bias
% dom       Named area or spherical cap semiopening angle
% Jmax      How many multitapers do you want to use for the spectral
%           estimation? [] for all (default)
% rotcoord  Center of spherical cap region if not northpole in degrees
%           [longitude colatitude], 0<=longitude<360, 0<=colatitude<=180
% Nspec     Global spectrum of the noise
%
% OUTPUT:
%
% spec      Local power spectrum for provided spherical-harmonic degrees L
% specvar 	Error bars for the local spectrum
%
% Last modified by plattner-at-alumni.ethz.ch, 09/15/2016

defval('Jmax',[])
defval('rotcoord',[])
defval('method',1)
defval('onorout',[])
defval('Nspec',[])
   
if ~ischar(data)   

Lwid=Ltap;%2*Ltap+1;    


% Evaluate the input coefficients on a grid
%[data,lon,lat]=plm2xyz(lmcosi);
data=data(:)';

% Get the multitaper coefficients
if ischar(dom) | isempty(rotcoord)
    if length(dom)==2
        [G,V]=glmalpharing(dom,Ltap);
    else
        [G,V]=glmalpha(dom,Ltap);
    end
else
    [G,V]=glmalphapto(dom,Ltap,rotcoord(1),rotcoord(2),[]);
end

% If you don't want to sum all Slepian functions,
% then you need to order them
if ~ischar(dom)
    [V,isrt]=sort(V,'descend');
    G=G(:,isrt);
end


% Define the number of Slepian functions
if isempty(Jmax)
    Jmax=(Ltap+1)^2;
else
    Jmax=min(Jmax,(Ltap+1)^2);            
end

% Evaluate the multitaper coefficients on the same grid as the data:
% First evaluate the spherical harmonics
% Remember: ylm doesn't have the (-1)^m phase, so just shift lon by
% 180
if Ltap==0
    Y=ylm(0,0,(90-lat)*pi/180,lon*pi/180+pi);   
    Y=Y(:)';
    G=full(G);
else
    Y=ylm([0 Ltap],[],(90-lat)*pi/180,lon*pi/180+pi);        
end
% evaluate the tapers
tapers=G'*Y;
gtimesd=tapers(1:Jmax,:).*repmat(data,Jmax,1);

% Get the spectrum for the sum of tapers*data
spec=zeros(Lmax+1,1);
sumV=sum(V(1:Jmax));                            

if ischar(dom)
    domarea=spharea(dom,0);
elseif length(dom)==1
    domarea=spharea(dom,1);
elseif length(dom)==2
    domarea=spharea(max(dom),1)-spharea(min(dom),1);
else
    error('Something is wrong with the domain')
end

%fact=4*pi;

fact=4*pi*domarea;

for alpha=1:Jmax
    % This is the product of the tapers times the data
    %gtimesd=tapers(alpha,:).*data;
    % This is the spectrum of each taper-data-product
    specprod=fact*...
        plm2spec(xyz2plm(...
        reshape(gtimesd(alpha,:),length(lat),length(lon)),Lmax)); 
    % If you want to see the individual tapered spectra: 
    if nargout>2
        spectap{alpha}=specprod;
    end
    % And normalize and sum them up
    spec=spec+V(alpha)/sumV*specprod;            
end        

if nargout<=2
    spectap=[];
end        

% Now get the error bars if you ask for them
if nargout>1            
    specvar=mtvar(spec,(0:Lmax)',Lwid,dom); 
else
    specvar=[];
end


% Now subtract the noise spectrum if provided
if ~isempty(Nspec)
    noiseLmax=length(Nspec)-1;
    M=mcouplings(Ltap,noiseLmax);
    disp('For the noise spectrum we just do eigenvalue weighted sum like in D.S. 2008 eq. (145)')
    spec=spec-M*Nspec(:);
end


varns={spec,specvar,spectap};
varargout=varns(1:nargout);

elseif strcmp(data,'demo1')
   disp('Demo 1: Comparison localspectrum to localspectrum2')
   
   Lmax=100;
   dom='namerica';
   Ltap=5;
   res=0.5;
   
   Jmax=3*round(spharea(dom)*(Ltap^2+1));
   
   lmcosi=plm2rnd(Lmax,rand(1));
   
      
   spec1=localspectrum(lmcosi,Ltap,dom,Jmax);
   

   [data,lon,lat]=plm2xyz(lmcosi,res);
   
   spec2=localspectrum2(data,lon,lat,Lmax,Ltap,dom,Jmax);
   
   plot(0:Lmax,spec1,'r-')
   hold on
   plot(0:Lmax,spec2,'--b')
   xlabel('degrees')
   ylabel('spectral power')
   legend('localspectrum','localspectrum2')
 
   %figure
   %plot(0:Lmax,spec1-spec2)
    
end
