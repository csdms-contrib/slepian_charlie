function [l,sdl,lmcosi,V,A]=mtap(lon,lat,data,L,TH,phi0,theta0,lmax,J)
colat=90-lat;
%L=round((N+1)*180/TH-1);
defval('L',8)
defval('J',(L+1)^2)
[Gar,V,n,J,phi0,theta0,omega,theta,phi]=...
    galphapto(TH,L,phi0,theta0,[],colat*pi/180,lon*pi/180,J);
length(V)
for i=1:length(V)
    disp(i)
    taper=reshape(Gar(i,:),length(lat),length(lon));
    tapdat=data.*taper;
    imagesc(tapdat);
    lmcosi(:,:,i)=xyz2plm(tapdat,lmax,[]);
    [sdl(:,i),l,bta,lfit,logy,logpm]=plm2spec(lmcosi(:,:,i),1);
end
A=2*pi*(1-cosd(TH));
