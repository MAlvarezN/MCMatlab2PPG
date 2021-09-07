% compute cummulative power
clear
close all
clc

h = figure ;
hold on

% fldResults = "..\Results\2021\08_24\bt2_t1\";
Nsteps = 100 ;
%% Collector information
load("E:\MCMatlab\Results\2021\09_02\model.mat")
% axis values
x = model.G.x ; % values used in the model
y = model.G.y ; % 
z = model.G.z ; % 
% Collector position [cm]
xcm = model.MC.LC.x ;
ycm = model.MC.LC.y ;
zcm = model.MC.LC.z ;
% indexes collector [pix]
[~,xi] = min(abs(x-xcm)) ;
[~,yi] = min(abs(y-ycm)) ;
[~,zi] = min(abs(z-zcm)) ;

% size collector 3.5mm side square?
side = 2.1 /10 ; % approx. 2.1 based on Figure 4 2020-Bonya

dx = model.G.dx ; % values used in the model
dy = model.G.dy ; % 
dz = model.G.dz ; % 

sdx = side / dx ; % pixels number for x dimension
sdy = side / dy ; % pixels number for y dimension
sdz = side / dz ; % pixels number for z dimension

% select cube
index_x = xi - floor( sdx / 2 ) : xi + floor( sdx / 2 ) ;
index_y = yi - floor( sdy / 2 ) : yi + floor( sdy / 2 ) ;
index_z = zi : zi + floor( sdz / 2 ) ;


%% trials
cummulativePowerAbsorTrial = zeros(5,Nsteps) ; 
for k = 1 : 1
    home
disp(["trial: " + k])
% Save current model: date + beam type
fldResults = "..\Results\2021\08_25\bt2_t" + k + "\";

cummlativePowerAbsor = zeros(1,Nsteps) ;
for j = 1 : Nsteps
%% load model    
    disp(["Step: " + j + " / " + Nsteps])
    
    nameModel = fldResults + "modelMC_pL1_" + j + "_" + Nsteps ;
    
    load(nameModel) % model structure

ymin = - .7 ;% model.G.Ly / 2 ; % supose that dont change
ymax = + .7 ;% model.G.Ly / 2 ;

ypos = linspace( ymin , ymax , Nsteps) ;

% MCorFMC = model.MC ;
mua_vec = [MCorFMC.mediaProperties.mua];

test2 = mua_vec(MCorFMC.M).*MCorFMC.NFR ;

% test2 = test2(index_x,index_z,index_z) ;

cum_ref = zeros( 1 , size(test2,2) ) ;
for i = 1:size(test2,2)
    cum_ref(i) = sum(sum( test2(:,i,:) )) ;
end

% cummulative
 cummlativePowerAbsor(j)  = sum(cum_ref) ;


    plot(ypos(1:j) , cummlativePowerAbsor(1:j) ,'-r')    
    drawnow
    xlabel('y[cm]')
    ylabel('total Power (sum(sum))')
    title('ppg?')
%     pause(.0001)
    
end

cummulativePowerAbsorTrial(k,:) = cummlativePowerAbsor ;

end
%% boxplot
figure,
boxplot(cummulativePowerAbsorTrial)

return
figure,
subplot(151),
    boxplot(cummulativePowerAbsorTrial(1:10,:))
subplot(152),
    boxplot(cummulativePowerAbsorTrial(1:20,:))
subplot(153),
    boxplot(cummulativePowerAbsorTrial(1:30,:))
subplot(154),
    boxplot(cummulativePowerAbsorTrial(1:40,:))
subplot(155),
    boxplot(cummulativePowerAbsorTrial(1:50,:))

figure,
subplot(151),
    plot(mean(cummulativePowerAbsorTrial(1:10,:)))
    axis square
subplot(152),
    plot(mean(cummulativePowerAbsorTrial(1:20,:)))
    axis square
subplot(153),
    plot(mean(cummulativePowerAbsorTrial(1:30,:)))
    axis square
subplot(154),
    plot(mean(cummulativePowerAbsorTrial(1:40,:)))
    axis square
subplot(155),
    plot(mean(cummulativePowerAbsorTrial(1:50,:)))
    axis square