% compute cummulative power
clear
close all
clc

h = figure ;
hold on

% fldResults = "..\Results\2021\08_24\bt2_t1\";
Nsteps = 100 ;

%% trials
cummulativePowerAbsorTrial = zeros(5,Nsteps) ; 
fldres = "09_11_700" ;
for k = 1 : 10
    home
disp(["trial: " + k])
% Save current model: date + beam type
fldResults = "..\Results\2021\"+fldres+"\bt2_t" + k + "\";

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
    ntitle = "ppg 700[nm] (10 trial)" ;
    title(ntitle)
%     pause(.0001)
    
end

cummulativePowerAbsorTrial(k,:) = cummlativePowerAbsor ;

end
%% boxplot
figure,
boxplot(cummulativePowerAbsorTrial)

load("D:\Michael_MCMatlab\Results\2021\"+fldres+"\model.mat")
figure,
% Distance (guess) for the Ellipses on geometry
    pL1 = .2 ;
    pL2 = -.3 ;

Nsteps = 100 ;
pL1values = linspace( - model.G.Ly/2 - (2* pL1 ) , model.G.Ly/2 + 2*( pL1 - pL2 ) , Nsteps ) ;
pL2values = pL1values - ( pL1 - pL2 ) ;

p = polyfit(pL1values,mean(cummulativePowerAbsorTrial),20) ;
    plot(mean(cummulativePowerAbsorTrial(:,1:i))) ; 
    hold on
    plot(polyval(p,pL1values(1:i)))
    title(ntitle)
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