%% Decription
% This example simulates a collimated top hat beam of radius 300 µm
% incident on skin, with some gel (water) on the top. This example is
% constructed identically to that on the mcxyz website, except that photons
% escape on all boundaries and the voxel grid is only 100x100x100:
% https://omlc.org/software/mc/mcxyz/
% 
% The found absorption distribution is then passed into the heat simulator,
% assuming the light is on for 5 pulses of 1 ms on time and 4 ms off time
% each, with 3 W of peak power. Some demonstration values of the Arrhenius
% E and A parameters for blood coagulation are used to calculate the
% distribution of coagulated blood. Temperature sensors outputs and movie
% generation is also demonstrated.

%

% close all; clear; clc

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 200; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 200; % Number of bins in the z direction
model.G.Lx                = 1.4; % [cm] x size of simulation cuboid
model.G.Ly                = 1.4; % [cm] y size of simulation cuboid
model.G.Lz                = 0.8; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc ; % Media properties defined as a function at the end of this file
model.G.geomFunc            = @geometryDefinition_BloodVessel ; 
                                                                % Function to use for defining the distribution of
                                                                % media in the cuboid. Defined at the end of this m file.

% ellipse parameters --------------------------------
pL1 = .2 ;
pL2 = -.3 ;
    
IDx  = .33 ; 
% pL1  = .2 ; 
pwd1 = .20  ; 

IDxd = .28 ; 
% pL2  = -.3 ; 
pwd2 = .40 ; 


model.G.geomFuncParams{1} = pL1 ; 
model.G.geomFuncParams{2} = pL2 ; 

model.G.geomFuncParams{3} = IDx ;
% model.G.geomFuncParams{4} = pL1 ;
model.G.geomFuncParams{5} = pwd1;

model.G.geomFuncParams{6} = IDxd;
% model.G.geomFuncParams{7} = pL2 ;
model.G.geomFuncParams{8} = pwd2;
% ---------------------------------------------------
% body mass parameters
bodymassindex = 25 ;
switch bodymassindex
    case 25
        vessel_diameter = 0.25 ;
        der_thick = 0.1 ; 
    case 30
        vessel_diameter = 0.275 ;
        der_thick = 0.1325 ; 
    case 35
        vessel_diameter = 0.30 ;
        der_thick = 0.175 ; 
    case 40
        vessel_diameter = 0.325 ;
        der_thick = 0.2125 ; 
    case 45
        vessel_diameter = 0.35 ;
        der_thick = 0.25 ; 
    otherwise
        disp('mmm... error!')
        return
end
model.G.geomFuncParams{9} = vessel_diameter ;
model.G.geomFuncParams{10} = der_thick ;

% plot(model,'G') ; % Dont need the plot in this part

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1 ; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true ; % Assumes all refractive indices are the same

model.MC.boundaryType             = 1; % 0: No escaping boundaries, 
                                       % 1: All cuboid boundaries are escaping, 
                                       % 2: Top cuboid boundary only is escaping
% lambda values: 660 - 750 - 800 - 850 - 900 - 950                                       
model.MC.wavelength               = 500 ; % [nm] Excitation wavelength,
                                         %  used for determination of optical properties for excitation 
                                         %  light
% model.G.mediaPropParams{1} = model.MC.wavelength ;

model.MC.beam.beamType            = 2; % 0: Pencil beam, 
                                       % 1: Isotropically emitting line or point source, 
                                       % 2: Infinite plane wave,
                                       % 3: Laguerre-Gaussian LG01 beam, 
                                       % 4: Radial-factorizable beam (e.g., a Gaussian beam), 
                                       % 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
                                       
model.MC.beam.NF.radialDistr      = 0; % Radial near field distribution - 
                                       %  0: Top-hat, 
                                       %  1: Gaussian, Array: Custom. Doesn't need to be normalized.
                                       
model.MC.beam.NF.radialWidth      = .03; % [cm] Radial near field 1/e^2 radius if top-hat or Gaussian 
                                         %  or half-width of the full distribution if custom
                                         
model.MC.beam.FF.radialDistr      = 0; % Radial far field distribution - 
                                       %  0: Top-hat, 
                                       %  1: Gaussian, 
                                       %  2: Cosine (Lambertian), 
                                       %  Array: Custom. Doesn't need to be normalized.
                                       
model.MC.beam.FF.radialWidth      = 0; % [rad] Radial far field 1/e^2 half-angle if top-hat or Gaussian 
                                       %  or half-angle of the full distribution if custom. 
                                       % For a diffraction limited Gaussian beam, this should be set to 
                                       %  model.MC.wavelength*1e-9/(pi*model.MC.beam.NF.radialWidth*1e-2))

ypos = 0.0 ;

model.MC.beam.xFocus              = .4 ; % [cm] x position of focus
model.MC.beam.yFocus              = ypos ; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

%% Light Collector
model.MC.useLightCollector      = true;
% model.MC.LC.x         = 0.1 ; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.x         = model.MC.beam.xFocus - 0.97 ;
model.MC.LC.y         = ypos ; % [cm] y position
model.MC.LC.z         = 0.0 ; % [cm] z position

model.MC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f         = .4;%.2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam      = .35; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize = .35; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA        = 0.866; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res       = 1; % X and Y resolution of light collector in pixels, only used for finite f

% model.MC.LC.tStart    = -1e-13; % [s] Start of the detection time-of-flight interval
% model.MC.LC.tEnd      = 5e-12; % [s] End of the detection time-of-flight interval
% model.MC.LC.nTimeBins = 30; % (Default: 0) Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved. 


%% Compute for each position


% Distance (guess) for the Ellipses on geometry
    pL1 = .2 ;
    pL2 = -.3 ;

Nsteps = 100 ;
pL1values = linspace( - model.G.Ly/2 - (2* pL1 ) , model.G.Ly/2 + 2*( pL1 - pL2 ) , Nsteps ) ;
pL2values = pL1values - ( pL1 - pL2 ) ;

% a = gcf ;
% figure(a.Number);

% p = polyfit(pL1values,mean(cummulativePowerAbsorTrial),20) ;

% start trials
for k = 1 : 10
rng(k)
    clc
disp(["trial: " + k])
% Save current model: date + beam type
fldResults0 = "..\Results\2021\09_11_500" ;
fldResults  = fldResults0 +         "\bt2_t" + k + "\";
mkdir(fldResults) % make the folder dont worry if the folder exists

time = zeros( 1 , Nsteps ) ;
for i = 1 : Nsteps
    disp(["Step: " + i + " / " + Nsteps ])
    model.G.geomFuncParams{1} = pL1values( i ) ; 
    model.G.geomFuncParams{2} = pL2values( i ) ; 

%     f=figure(1);
%     f.Position = [1   100   895   648];
%         plot(model,'G') ; 
%         drawnow
%     g=figure(2);
%     g.Position = [900   100   711   648];
%         plot(mean(cummulativePowerAbsorTrial(:,1:i))) ; 
%         drawnow
%         hold on
%         plot(polyval(p,pL1values(1:i)))
%         title('PPG')
% % 
tic    
% Execution, do not modify the next line:
model = runMonteCarlo(model);
plot(model,'MC');
time( i ) = toc ;

if k==1 && i==1
    disp(" Save one model to save all parameters...")
    nameModel = fldResults0 + "\model" ;
    save(nameModel,'model') ;
end
    

% Save current model
nameModel = fldResults + "modelMC_pL1_" + i + "_" + Nsteps ;
MCorFMC = model.MC ;
save(nameModel,'MCorFMC') ;

end

% Save time execution
nameTimeModel = fldResults + "time_modelMC_pL1" ;
save(nameTimeModel,'time') ;

end % trials

% reflectance_vessel

%% Post-processing

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of Ginput. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = geometryDefinition_BloodVessel(X,Y,Z,parameters)
% Blood vessel example:
zsurf     = 0.005; % cm
epd_thick = 0.01;
% der_thick = 0.1 ; 
vessel_thick = 0.02 ; 
% vessel_diameter = 0.25 ;

try
    vessel_diameter =  parameters{9} ;
    der_thick       =  parameters{10} ; 
catch
    vessel_diameter = 0.25 ;
    der_thick = 0.1 ; 
end
vesselradiusOut  = vessel_diameter / 2 ;
vesselradius  = vesselradiusOut - vessel_thick ;
vesseldepth = 0.4;

M = ones(size(X)); % fill background with water (gel)
M(Z > zsurf) = 6; % rEpidermis
M(Z > zsurf + epd_thick) = 2; % dermis
M(Z > zsurf + epd_thick + der_thick) = 5 ; % subcutis
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradiusOut^2) = 4 ; % vessel
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 3 ; % blood

% ellipsoid eqtuations ----------------------------------------------------

try
    pL1 = parameters{1} ;
    pL2 = parameters{2} ;
    

    IDx  = parameters{3} ;
%     pL1  = parameters{4} ;
    pwd1 = parameters{5} ;

    IDxd = parameters{6} ;
%     pL2  = parameters{7} ;
    pwd2 = parameters{8} ;

catch
    pL1 = .2 ;
    pL2 = -.3 ;
    
    IDx  = .33 ; 
%     pL1  = .2 ; 
    pwd1 = .33 ; 

    IDxd = .28 ; 
%     pL2  = -.3 ; 
    pwd2 = .28 ; 

end
d    = zsurf + vesseldepth ;

cond_Eq2 = (X./(IDx/2)).^2 + ((Z-d)./(IDx/2)).^2 + ((Y+pL1)./(pwd1)).^2 <= 1 ;
M( cond_Eq2 ) = 4 ; % vessel

cond_Eq3 = (X./(IDxd/2)).^2 + ((Z-d)./(IDxd/2)).^2 + ((Y+pL2)./(pwd2)).^2 <= 1 ;
M( cond_Eq3 ) = 4 ; % vessel

% Adjust the blood in vessel
cond_Eq22 = (X./(IDx/2 - vessel_thick)).^2 + ((Z-d)./(IDx/2 - vessel_thick)).^2 + ((Y+pL1)./(pwd1 - vessel_thick)).^2 <= 1 ;
M( cond_Eq22 ) = 3 ; % blood
cond_Eq32 = (X./(IDxd/2 - vessel_thick)).^2 + ((Z-d)./(IDxd/2 - vessel_thick)).^2 + ((Y+pL2)./(pwd2 - vessel_thick)).^2 <= 1 ;
M( cond_Eq32 ) = 3 ; % blood

% add blood
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 3 ; % blood


end

%% Media Properties function
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
load spectralLIB.mat
MU(:,1) = interp1(nmLIB,muaoxy,wavelength);
MU(:,2) = interp1(nmLIB,muadeoxy,wavelength);
MU(:,3) = interp1(nmLIB,muawater,wavelength);
MU(:,4) = interp1(nmLIB,muamel,wavelength);

j=1;
mediaProperties(j).name  = 'water';
mediaProperties(j).mua   = 0.00036;
mediaProperties(j).mus   = 10;
mediaProperties(j).g     = 1.0;
mediaProperties(j).n     = 1.3;
mediaProperties(j).VHC   = 4.19;
mediaProperties(j).TC    = 5.8e-3;

% j=2;
% mediaProperties(j).name  = 'epidermis';
% B = 0;
% S = 0.75;
% W = 0.75;
% Me = 0.03;
% musp500 = 40;
% fray    = 0.0;
% bmie    = 1.0;
% gg      = 0.90;
% musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
% X = [B*S B*(1-S) W Me]';
% mediaProperties(j).mua = MU*X;
% mediaProperties(j).mus = musp/(1-gg);
% mediaProperties(j).g   = gg;
% mediaProperties(j).n   = 1.3;
% mediaProperties(j).VHC = 3391*1.109e-3;
% mediaProperties(j).TC  = 0.37e-2;

% j=3;
% mediaProperties(j).name = 'dermis';
% B = 0.002;
% S = 0.67;
% W = 0.65;
% Me = 0;
% musp500 = 42.4;
% fray    = 0.62;
% bmie    = 1.0;
% gg      = 0.90;
% musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
% X = [B*S B*(1-S) W Me]';
% mediaProperties(j).mua = MU*X;
% mediaProperties(j).mus = musp/(1-gg);
% mediaProperties(j).g   = gg;
% mediaProperties(j).n   = 1.3;
% mediaProperties(j).VHC = 3391*1.109e-3;
% mediaProperties(j).TC  = 0.37e-2;

% j=4;
% mediaProperties(j).name  = 'blood';
% B       = 1.00;
% S       = 0.75;
% W       = 0.95;
% Me      = 0;
% musp500 = 10;
% fray    = 0.0;
% bmie    = 1.0;
% gg      = 0.90;
% musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
% X = [B*S B*(1-S) W Me]';
% mediaProperties(j).mua = MU*X;
% mediaProperties(j).mus = musp/(1-gg);
% mediaProperties(j).g   = gg;
% mediaProperties(j).n   = 1.3;
% mediaProperties(j).VHC = 3617*1.050e-3;
% mediaProperties(j).TC  = 0.52e-2;
% mediaProperties(j).E   = 422.5e3; % J/mol    PLACEHOLDER DATA ONLY
% mediaProperties(j).A   = 7.6e66; % 1/s        PLACEHOLDER DATA ONLY


switch wavelength
    
% 660 nm ------------------------------------------------------------------
    case 660

j=2;
mediaProperties(j).name = 'dermis';
mediaProperties(j).mua = 0.5453 ;
mediaProperties(j).mus = 208.6 ;
mediaProperties(j).g   = 0.7 ;
mediaProperties(j).n   = 1.47 ;
% mediaProperties(j).VHC = 3391*1.109e-3;
% mediaProperties(j).TC  = 0.37e-2;

j=3;
mediaProperties(j).name  = 'blood';
mediaProperties(j).mua = 2.026 ;
mediaProperties(j).mus = 75.76 ;
mediaProperties(j).g   = 0.9 ;
mediaProperties(j).n   = 1.4 ;
% mediaProperties(j).VHC = 3617*1.050e-3;
% mediaProperties(j).TC  = 0.52e-2;
% mediaProperties(j).E   = 422.5e3; % J/mol    PLACEHOLDER DATA ONLY
% mediaProperties(j).A   = 7.6e66; % 1/s        PLACEHOLDER DATA ONLY

j = 5 ;
mediaProperties(j).name  = 'subcutis';
mediaProperties(j).mua   = 0.0001 ;
mediaProperties(j).mus   = 249.7 ;
mediaProperties(j).g     = 0.7 ;
mediaProperties(j).n     = 1.47 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

j = 4 ;
mediaProperties(j).name  = 'vessel';
mediaProperties(j).mua   = 0.8 ;
mediaProperties(j).mus   = 230 ;
mediaProperties(j).g     = 0.9 ;
mediaProperties(j).n     = 1.4 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

j = 6 ;
mediaProperties(j).name  = 'rEpidermis';
mediaProperties(j).mua   = 0.3442 ;
mediaProperties(j).mus   = 121.2 ;
mediaProperties(j).g     = 0.7 ;
mediaProperties(j).n     = 1.47 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

% 890 nm ------------------------------------------------------------------
    case 890

j=2;
mediaProperties(j).name = 'dermis';
mediaProperties(j).mua = 0.2459 ;
mediaProperties(j).mus = 116.7 ;
mediaProperties(j).g   = 0.7 ;
mediaProperties(j).n   = 1.47 ;
% mediaProperties(j).VHC = 3391*1.109e-3;
% mediaProperties(j).TC  = 0.37e-2;

j=3;
mediaProperties(j).name  = 'blood';
mediaProperties(j).mua =  6.32 ;
mediaProperties(j).mus = 56.18 ;
mediaProperties(j).g   =  0.9 ;
mediaProperties(j).n   =  1.4 ;
% mediaProperties(j).VHC = 3617*1.050e-3;
% mediaProperties(j).TC  = 0.52e-2;
% mediaProperties(j).E   = 422.5e3; % J/mol    PLACEHOLDER DATA ONLY
% mediaProperties(j).A   = 7.6e66; % 1/s        PLACEHOLDER DATA ONLY

j = 5 ;
mediaProperties(j).name  = 'subcutis';
mediaProperties(j).mua   = 0.0217 ;
mediaProperties(j).mus   = 189.8 ;
mediaProperties(j).g     = 0.7 ;
mediaProperties(j).n     = 1.47 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

j = 4 ;
mediaProperties(j).name  = 'vessel';
mediaProperties(j).mua   = 0.8 ;
mediaProperties(j).mus   = 230 ;
mediaProperties(j).g     = 0.9 ;
mediaProperties(j).n     = 1.4 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

j = 6 ;
mediaProperties(j).name  = 'rEpidermis';
mediaProperties(j).mua   = 0.3184 ;
mediaProperties(j).mus   = 224.7 ;
mediaProperties(j).g     = 0.7 ;
mediaProperties(j).n     = 1.47 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

% 500 nm ------------------------------------------------------------------
    case 500

j=2;
mediaProperties(j).name = 'dermis';
mediaProperties(j).mua = 1.51 ;
mediaProperties(j).mus = 424 ;
mediaProperties(j).g   = 0.9 ;
mediaProperties(j).n   = 1.4 ;
j=3;
mediaProperties(j).name  = 'blood';
mediaProperties(j).mua = 1.12 ;%112 ?
mediaProperties(j).mus = 100 ;
mediaProperties(j).g   = 0.9 ;
mediaProperties(j).n   = 1.4 ;
j = 5 ;
mediaProperties(j).name  = 'subcutis';
mediaProperties(j).mua   = 0.0001 ;
mediaProperties(j).mus   = 350 ;
mediaProperties(j).g     = 0.9 ;
mediaProperties(j).n     = 1.4 ;
j = 4 ;
mediaProperties(j).name  = 'vessel';
mediaProperties(j).mua   = 0.8 ;
mediaProperties(j).mus   = 230 ;
mediaProperties(j).g     = 0.9 ;
mediaProperties(j).n     = 1.4 ;
j = 6 ;
mediaProperties(j).name  = 'rEpidermis';
mediaProperties(j).mua   = 2.16 ;
mediaProperties(j).mus   = 400 ;
mediaProperties(j).g     = 0.9 ;
mediaProperties(j).n     = 1.4 ;
% -------------------------------------------------------------------------
    otherwise
        disp('mmm...')
        disp('Unknown parameters for selected wavelength')
        return
end


end
