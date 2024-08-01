function varargout=localspectrum(lmcosi,Ltap,dom,Jmax,rotcoord,Nspec,method,optn,rplanet)
% [spec,spectap,V,sig,specvar]=localspectrum(lmcosi,Ltap,dom,Jmax,rotcoord,Nspec,method,optn,rplanet)
%
% Calculates the local multitaper spectrum using the formula described by
% Dahlen & Simons (2008), eq. 130
%
% INPUT:
%
% lmcosi    Spherical harmonic coefficients of the function for which you
%           want to calculate the local spectrum
% Ltap      Maximum spherical-harmonic degree for the multitapers with
%           which you want to calculate the local spectrum.
% dom       Named area or spherical cap semiopening angle
% Jmax      How many multitapers do you want to use for the spectral
%           estimation? [] for all (default)
% rotcoord  Center of spherical cap region if not northpole in degrees
%           [longitude colatitude], 0<=longitude<360, 0<=colatitude<=180
% Nspec     Global spectrum of the noise
% method 	  There are different ways to do this:
% 			    1 evaluate, multiply with tapers, get spectrum, sum
% 		 	    2 Use analytical expression for integral of three spherical-harmonics
% 			    3 Appendix of Dahlen & Simons (2008)
% 			    4 plm2spec(GJ*GJ'c)
% optn      0 (default) use division by (2*l+1) spectral normalization and 
%             divide by area*4*pi
%           1 use multiplication by (l+1) spectral normalization, no
%             division
%           2 normalize to obtain Mauersberger-Lowes spectrum
% rplanet   Planet radius. Only required for Mauersberger-Lowes spectrum.  
%  
% OUTPUT:
%
% spec        Local power spectrum for provided spherical-harmonic degrees L
% spectap     The individual tapered spectra
% V           concentration values
% sig         standard deviation from the individual single taper spectra
% specvar 	  Error bars for the local spectrum (currently broken)
%
% Last modified by plattner-at-alumni.ethz.ch, 05/01/2024

defval('Jmax',[])
defval('rotcoord',[])
defval('method',1)
defval('onorout',[])
defval('Nspec',[])
defval('optn',0)

if ~ischar(lmcosi)    
    
Lwid=Ltap;%2*Ltap+1;    

if optn==1
    fact=4*pi;
    specnorm=1;
else
    if ischar(dom)
        domarea=spharea(dom,0);
    elseif length(dom)==1
        domarea=spharea(dom,1);
    elseif length(dom)==2
        domarea=spharea(max(dom),1)-spharea(min(dom),1);
    else
        error('Something is wrong with the domain')
    end
    
    fact=4*pi*domarea;
    specnorm=2;
end

switch method
	case 1
        % Turn coefs into lmcosi       
        Lmax=lmcosi(end,1);
        
        % Evaluate the input coefficients on a grid
        [data,lon,lat]=plm2xyz(lmcosi);
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
        % tapers=G'*Y;
        % gtimesd=tapers(1:Jmax,:).*repmat(data,Jmax,1);
        % The following should be faster and save memory
        GJ=G(:,1:Jmax);
        clear G
        tapers=GJ'*Y;
        gtimesd=tapers.*repmat(data,Jmax,1);
        
        % Get the spectrum for the sum of tapers*data
        spec=zeros(Lmax+1,1);
        sumV=sum(V(1:Jmax));                            
        
     
        
        for alpha=1:Jmax
            % This is the product of the tapers times the data
            %gtimesd=tapers(alpha,:).*data;
            % This is the spectrum of each taper-data-product
            specprod=fact*...
                plm2spec(xyz2plm(...
                reshape(gtimesd(alpha,:),length(lat),length(lon)),Lmax),specnorm); 
            % If you want to see the individual tapered spectra: 
            if nargout>1
                spectap{alpha}=specprod;
            end
            % And normalize and sum them up
            spec=spec+V(alpha)/sumV*specprod;            
        end        
           
        if nargout<=1
            spectap=[];
            specvar=[];
            sig = [];        
        else
            warning('Since April 8, 2024: Second output is cell of individually tapered spec and 4th output is error bars')
        end        
        
        % Now get the error bars if you ask for them
        if nargout>4            
            try
                specvar=mtvar(spec,(0:Lmax)',Lwid,dom); 
                specvar=specvar(:);
            catch
                keyboard
                %warning('problems with wignercycle, thus mtvar is currently out of order')
                %specvar = [];

                %warning('calculate first wigner symbols. This may take a while but only the first time')
                %wignercycle(2*Lmax);
                %specvar=mtvar(spec,(0:Lmax)',Lwid,dom);
                %keyboard
            end
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
        
	otherwise
		error('Not yet implemented')
end

if optn==2
  %% Normalize for Mauersberger-Lowes spectrum
  Lmax=max(lmcosi(:,1));
  ls=(0:Lmax)';
  spec=spec.*(ls+1).*(2*ls+1).^2/rplanet^2;

  if nargout>1
    for alpha=1:Jmax
        spectap{alpha}=spectap{alpha}.*(ls+1).*(2*ls+1).^2/rplanet^2;
    end
  end
end

if nargout>1
    SV = cell2mat(spectap);
    sig = std(SV,abs(V(1:Jmax)),2); % For some regions, the eigenvalues can be around -eps.
end

varns={spec,spectap,V(:),sig,specvar};
varargout=varns(1:nargout);

elseif strcmp(lmcosi,'demo1')
    % Check if it does what we want   
    Lmax=100;    
    dom=90;%'africa';
    Jmax=[];
    rotcoord=[];
    beta=-0;
    spectap=[];
    if ~ischar(dom)
        Ltap=round(2*180/dom);
    end
    %Ltap=1;
    
    % Make a field that globally has a spectrum
    lmcosi=plm2rnd(Lmax,beta);
      
    % Get the multitaper local spectrum estimation
    [spec,specvar,spectap]=localspectrum(lmcosi,Ltap,dom,Jmax,rotcoord);   
    
    fact=1;%(4*pi)
    
    if isempty(Jmax)
        Jmax=(Ltap+1)^2;
    end
    
    % This is the true global spectrum:
    truespec=plm2spec(lmcosi)/fact;      
    
    % Plot the result
    lw=1.5;
    lw2=0.2;
    plot(0:Lmax,truespec,'k','LineWidth',lw)
    hold on
    plot(0:(Lmax-Ltap),spec(1:end-Ltap),'b--','LineWidth',lw)
    %uncor=round(Ltap/2):2*Ltap:(Lmax-Ltap);
    uncor=[round(Ltap/2),2*Ltap:2*Ltap:Lmax-Ltap]; % See mcoupling
    %plot(uncor,spec(uncor+1),'bo')
    if ~isempty(spectap)
    for alpha=1:Jmax
        plot(0:(Lmax-Ltap),spectap{alpha}(1:end-Ltap),'-','LineWidth',lw2,'color',[0.7 0.7 0.7])
    end
    end
    errorbar(uncor,spec(uncor+1),sqrt(specvar(uncor+1)),'bo','LineWidth',lw)    
    legend('true','multitaper','single taper')
    xlabel('L')
    ylabel('power spectrum')
    
    % And again
    plot(0:Lmax,truespec,'k','LineWidth',lw)
    plot(0:(Lmax-Ltap),spec(1:end-Ltap),'b--','LineWidth',lw)    
    
    ylim([0 max(truespec)*1.2])
    
    title(sprintf('Estimation for L=%d and Jmax=%d',Ltap,Jmax))
    hold off
    
    if ischar(dom)
        print(sprintf('%s_Ltap%d',dom,Ltap),'-dpdf')
    else
        print(sprintf('TH%d_Ltap%d',dom,Ltap),'-dpdf')
    end
    
    disp('Also show the espected value given the truth: M*truespec')
    if beta==0      
        ylim([0 2/fact])
    end
    
elseif strcmp(lmcosi,'demo2')
    % Try check if error bars are ok if running many realizations of this estimation    
     
    Lmax=100;   
    dom=179;%'africa';
    Jmax=[];
    rotcoord=[];    
    ntry=10;          
    
    % To get a good Ltap for spherical caps:
    if ~ischar(dom)
        Ltap=round(180/dom);
        %Lwid=round(2*180/dom);
    end
    %Ltap=1;
    
    Lwid=Ltap;%2*Ltap+1;
    
    % The true spectrum from which we randomly pick:
    beta=0;
    truespec=((1:Lmax+1).^beta)';
    
    % Plot it      
    lw2=1.5;
    plot(0:Lmax,truespec,'k','LineWidth',lw2)  
    hold on
    uncor=[round(Lwid/2),2*Lwid:2*Lwid:Lmax-Lwid]; % See mcoupling
    % Get the error bars from the true spectrum
    %specvar=mtvar(truespec,(0:Lmax)',Ltap,dom);
    specvar=mtvar(truespec,(0:Lmax)',Lwid,dom);
    % Plot the error bars
    errorbar(uncor,truespec(uncor+1),2*sqrt(specvar(uncor+1)),...
        'o','LineWidth',lw2,'color',[0.3 0.3 0.3])
    errorbar(uncor,truespec(uncor+1),sqrt(specvar(uncor+1)),...
        'ko','LineWidth',lw2)    
    
    meanspec=zeros(size(truespec));
    for i=1:ntry
        % Make a random field that globally has a spectrum
        lmcosi=spec2rnd(Lmax,truespec); 
        % Get the multitaper local spectrum estimation
        spec=localspectrum(lmcosi,Ltap,dom,Jmax,rotcoord);   
        meanspec=meanspec+spec;
        plot(0:(Lmax-Ltap),spec(1:end-Ltap),'-','LineWidth',lw2,'color',[0.7 0.7 0.7])
        plot(uncor,spec(uncor+1),'o','color',[0.7 0.7 0.7])
    end
    meanspec=meanspec/ntry;

    % And again to make it visible
    plot(0:Lmax,truespec,'k','LineWidth',lw2)   
    errorbar(uncor,truespec(uncor+1),sqrt(specvar(uncor+1)),...
        'ko','LineWidth',lw2)
    plot(0:Lmax,meanspec,'b','LineWidth',lw2)  
    
    if beta==0
        ylim([0 2])
    end
    hold off
    
    title(sprintf('TH=%d, Ltap=%d',dom,Ltap))
    xlabel('l')
    ylabel('S_l')
    
end

