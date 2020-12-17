function[E] = fretsimulation(extCoeFA,n,QYD,ltD,h,r,D,d1pos,d2pos)
h = input('plz enter fiber height');
% cyllinder height 
r= input ('plz enter fiber radius');
%cyllinder radius
[X,Y,Z] = cylinder(r);
Z = Z*h;
%build cyllinder
d1pos=input('plz enter [theta1 rho1 z1] for dipole1 ');
d2pos=input('plz enter [theta2 rho2 z2] for dipole1');
[x1,y1,z1] = pol2cart(d1pos(1,1),d1pos(1,2),d1pos(1,3));
[x2,y2,z2] = pol2cart(d2pos(1,1),d2pos(1,2),d2pos(1,3));
D=sqrt(abs((x1-x2)^2)+abs((y1-y2)^2)+abs((z1-z2)^2)); 
%distance between dipoles
N_a = 6.022e23;%avogado number
n   = input('plz insert refrctive index'); %refIndex
  % contant for Forster integral calculation
f1 = 9000 * log(10) / (128 * pi^5 * N_a); % [mol cm^3]
QYD =input('plz insert Quantum yield');
teta_d= input('please insert tetaD--->');
%tetaD is the angle between donnor vector and vector tht join donnor and
%acceptor,angle must be in degree
teta_a= input('please inpute tetaA--->');
%tetaA is the angle between acceptor vector and vector tht join donnor and
%acceptor,angle must be in degree
teta_t=input('please insert tetaT--->');
%tetaT is the angle between emission transition dipole of donnor and
%transition absorption of the acceptor,angle must be in degree
k2=(cos(teta_t)-3*cos(teta_d)*cos(teta_a))^2;%orientation factor
ltD=input('please enter life time of donnor');%lifetime of donor
load('DY490PBSSpecData (3).mat')
% abs spec of the acceptor 
absA = specdata(:,2); 
% emission spec of the donoer
emD  = specdata(:,3); 
% common wavelength axis
wavelength = specdata(:,1); % [nm] 
  %   wavelength: wavelength axis
assert(iscolumn(absA),'Absorption of the acceptor must be a colum vector')
assert(iscolumn(emD),'Emission of the donor must be a colum vector')
assert(iscolumn(wavelength),'Wavelength axis must be a colum vector')
assert(and(length(absA)==length(wavelength),length(emD)==length(wavelength)),...
'abs and emission spec must be on the same wavelength axis');
  % delta in the wavelength axis
norAbsA = absA ./ max(absA);
 % epsilon from Forster equation
extCoeFA=input('please enter extcoefA');
epsilon      = norAbsA .*extCoeFA; % [mol^-1 cm^-1]
  % sum normalized fluorescence spectra
fluoNorm    = (emD.*wavelength) ./ (sum(emD.*wavelength)); %[dimensionless]
 %   Spectral overlap parameter - J in Forster theory
J   = sum(fluoNorm .* epsilon .* (wavelength.^4)); % [mol^-1 cm^-1 nm^4]
f2 = 1e14; % [nm^2 * cm^-2] This factor takes into account the units of J
Rf = (f1 * f2 * ((QYD * J * k2) ./ (n^4)))^(1/6);% [nm];
E=((Rf).^6/(Rf.^6+D.^6));
figure,scatter(E,D);
title('foster energy');
xlabel('dipole distance') ;
ylabel('foster energy') ;
assignin('base','E',E);
assignin('base','Rf',Rf);
assignin('base','k2',k2);
assignin('base','D',D);
figure,mesh(X,Y,Z);
hold on
scatter3(d1pos(1,1),d1pos(1,2),d1pos(1,3));
hold on
scatter3(d2pos(1,1),d2pos(1,2),d2pos(1,3));
end
