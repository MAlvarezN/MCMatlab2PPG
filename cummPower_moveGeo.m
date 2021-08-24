%% compute cummulative power

h = figure ;

fldResults = "..\Results\2021\08_24\";
Nsteps = 100 ;
cummlativePowerAbsor = zeros(1,Nsteps) ;
for j = 1 : Nsteps
    disp(["Step: " + j + " / " + Nsteps])
    
    nameModel = fldResults + "modelMC_pL1_" + j + "_" + Nsteps ;
    
    load(nameModel) % model structure

ymin = - model.G.Ly / 2 ; % supose that dont change
ymax = + model.G.Ly / 2 ;

ypos = linspace( ymin , ymax , Nsteps) ;

MCorFMC = model.MC ;
mua_vec = [MCorFMC.mediaProperties.mua];

test2 = mua_vec(MCorFMC.M).*MCorFMC.NFR ;

cum_ref = zeros( 1 , size(test2,2) ) ;
for i = 1:size(test2,2)
    cum_ref(i) = sum(sum( test2(:,i,:) )) ;
end

% cummulative
 cummlativePowerAbsor(j)  = sum(cum_ref) ;


    plot(ypos(1:j) , cummlativePowerAbsor(1:j) )    
    drawnow
    xlabel('y[cm]')
    ylabel('total Power (sum(sum))')
    title('ppg?')
    pause(.0001)
    
end
%     
%   
% 
% ypos = linspace( ymin , ymax , 89) ;
% cummlativePowerAbsor = zeros(size(ypos)) ;
% figure,
% for j = 1 : length(ypos)
% 
%     model.MC.beam.yFocus = ypos(j) ; % [cm] y position of focus
%     model.MC.LC.y        = ypos(j) ; % [cm] y position
% 
%     % Execution, do not modify the next line:
%     model = runMonteCarlo(model);
% 
%     % plot(model,'MC');
% 
% 
%     reflectance_vessel
%     cummlativePowerAbsor(j)  = sum(cum_ref) ;
% 
%     plot(ypos(1:j) , cummlativePowerAbsor(1:j) )
%     drawnow
%     pause(.0001)
% end
