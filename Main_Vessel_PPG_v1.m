%% Decription
% This example simulates a collimated top hat beam of radius 300 �m
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
close all; clear; clc

%% Geometry definition
model = MCmatlab.model;

model.G.nx                = 200; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 200; % Number of bins in the z direction
model.G.Lx                = 1.4; % [cm] x size of simulation cuboid
model.G.Ly                = 1.4; % [cm] y size of simulation cuboid
model.G.Lz                = 0.8; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition_BloodVessel; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

plot(model,'G');
% return
%% Monte Carlo simulation
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same

model.MC.boundaryType             = 1; % 0: No escaping boundaries, 
                                       % 1: All cuboid boundaries are escaping, 
                                       % 2: Top cuboid boundary only is escaping
                                       
model.MC.wavelength               = 660; % [nm] Excitation wavelength,
                                         %  used for determination of optical properties for excitation light

model.MC.beam.beamType            = 4; % 0: Pencil beam, 
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

model.MC.beam.xFocus              = 0; % [cm] x position of focus
model.MC.beam.yFocus              = 0; % [cm] y position of focus
model.MC.beam.zFocus              = 0; % [cm] z position of focus
model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis

%% Light Collector
model.MC.useLightCollector      = true;
model.MC.LC.x         = 0.1 ; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y         = 0.1 ; % [cm] y position
model.MC.LC.z         = 0.0 ; % [cm] z position

model.MC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f         = .2;%.2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res       = 1;%50; % X and Y resolution of light collector in pixels, only used for finite f

% model.MC.LC.tStart    = -1e-13; % [s] Start of the detection time-of-flight interval
% model.MC.LC.tEnd      = 5e-12; % [s] End of the detection time-of-flight interval
% model.MC.LC.nTimeBins = 30; % (Default: 0) Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved. 


% Execution, do not modify the next line:
model = runMonteCarlo(model);

plot(model,'MC');

%% Heat simulation
% model.MC.P                   = 3; % [W] Incident pulse peak power (in case of infinite plane waves, only the power incident upon the cuboid's top surface)
% 
% model.HS.useAllCPUs          = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% model.HS.makeMovie           = true; % Requires silentMode = false.
% model.HS.largeTimeSteps      = false; % (Default: false) If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.
% 
% model.HS.heatBoundaryType    = 0; % 0: Insulating boundaries, 1: Constant-temperature boundaries (heat-sinked)
% model.HS.durationOn          = 0.001; % [s] Pulse on-duration
% model.HS.durationOff         = 0.004; % [s] Pulse off-duration
% model.HS.durationEnd         = 0.02; % [s] Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train
% model.HS.Tinitial            = 37; % [deg C] Initial temperature
% 
% model.HS.nPulses             = 5; % Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use nPulses = 1.
% 
% model.HS.plotTempLimits      = [37 100]; % [deg C] Expected range of temperatures, used only for setting the color scale in the plot
% model.HS.nUpdates            = 100; % Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)
% model.HS.slicePositions      = [.5 0.6 1]; % Relative slice positions [x y z] for the 3D plots on a scale from 0 to 1
% model.HS.tempSensorPositions = [0 0 0.038
%                                 0 0 0.04
%                                 0 0 0.042
%                                 0 0 0.044]; % Each row is a temperature sensor's absolute [x y z] coordinates. Leave the matrix empty ([]) to disable temperature sensors.
% % 
% % % Execution, do not modify the next line:
% % % model = simulateHeatDistribution(model);
% % % 
% % % plot(model,'HS');

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
der_thick = 0.1 ; 
vessel_thick = 0.02 ; 
vessel_diameter = 0.25 ;
vesselradiusOut  = vessel_diameter / 2 ;
vesselradius  = vesselradiusOut - vessel_thick ;
vesseldepth = 0.4;
M = ones(size(X)); % fill background with water (gel)
M(Z > zsurf) = 7; % rEpidermis
M(Z > zsurf + epd_thick) = 3; % dermis
M(Z > zsurf + epd_thick + der_thick) = 5 ; % subcutis
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradiusOut^2) = 6 ; % vessel
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 4 ; % blood

% ellipsoid eqtuations ----------------------------------------------------
d    = zsurf + vesseldepth ; 

IDx  = .33 ; 
pL1  = .2 ; 
pwd1 = .33 ; 

IDxd = .28 ; 
pL2  = -.3 ; 
pwd2 = .28 ; 

cond_Eq2 = (X./(IDx/2)).^2 + ((Z-d)./(IDx/2)).^2 + ((Y+pL1)./(pwd1)).^2 <= 1 ;
M( cond_Eq2 ) = 6 ; % vessel

cond_Eq3 = (X./(IDxd/2)).^2 + ((Z-d)./(IDxd/2)).^2 + ((Y+pL2)./(pwd2)).^2 <= 1 ;
M( cond_Eq3 ) = 6 ; % vessel

% Adjust the blood in vessel
cond_Eq22 = (X./(IDx/2 - vessel_thick)).^2 + ((Z-d)./(IDx/2 - vessel_thick)).^2 + ((Y+pL1)./(pwd1 - vessel_thick)).^2 <= 1 ;
M( cond_Eq22 ) = 4 ; % blood
cond_Eq32 = (X./(IDxd/2 - vessel_thick)).^2 + ((Z-d)./(IDxd/2 - vessel_thick)).^2 + ((Y+pL2)./(pwd2 - vessel_thick)).^2 <= 1 ;
M( cond_Eq32 ) = 4 ; % blood

% add blood
M(X.^2 + (Z - (zsurf + vesseldepth)).^2 < vesselradius^2) = 4 ; % blood


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

j=2;
mediaProperties(j).name  = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
Me = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X;
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
mediaProperties(j).VHC = 3391*1.109e-3;
mediaProperties(j).TC  = 0.37e-2;

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

% 660 nm ------------------------------------------------------------------

j=3;
mediaProperties(j).name = 'dermis';
mediaProperties(j).mua = 0.5453 ;
mediaProperties(j).mus = 208.6 ;
mediaProperties(j).g   = 0.7 ;
mediaProperties(j).n   = 1.47 ;
% mediaProperties(j).VHC = 3391*1.109e-3;
% mediaProperties(j).TC  = 0.37e-2;

j=4;
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

j = 6 ;
mediaProperties(j).name  = 'vessel';
mediaProperties(j).mua   = 0.8 ;
mediaProperties(j).mus   = 230 ;
mediaProperties(j).g     = 0.9 ;
mediaProperties(j).n     = 1.4 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;

j = 7 ;
mediaProperties(j).name  = 'rEpidermis';
mediaProperties(j).mua   = 0.3442 ;
mediaProperties(j).mus   = 121.2 ;
mediaProperties(j).g     = 0.7 ;
mediaProperties(j).n     = 1.47 ;
% mediaProperties(j).VHC   = 4.19;
% mediaProperties(j).TC    = 5.8e-3;




end