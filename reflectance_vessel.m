
MCorFMC = model.MC ;
mua_vec = [MCorFMC.mediaProperties.mua];

test2 = mua_vec(MCorFMC.M).*MCorFMC.NFR ;

cum_ref = zeros( 1 , size(test2,2) ) ;
for i = 1:size(test2,2)
    cum_ref(i) = sum(sum( test2(:,i,:) )) ;
end

% figure,
%     plot(cum_ref)
%     xlabel('y-axis []')
%     ylabel('Cummulative Power Absorption []')

sum(cum_ref)
return
 %% Make power absorption plot
  % Calculate 3D absorption distribution, which may be FR or T dependent
  mua_vec = [MCorFMC.mediaProperties.mua];
  h_f = plotVolumetric.plotVolumetric(3 + figNumOffset,G.x,G.y,G.z,mua_vec(MCorFMC.M).*MCorFMC.NFR,'MCmatlab_fromZero');
  set(h_f,'WindowStyle','Docked');
  h_f.Name = ['Normalized ' fluorescenceOrNothing 'absorption'];
  title(['Normalized ' fluorescenceOrNothing 'absorbed power per unit volume [W/cm^3/W.incident]'])
  
  
  [0.000360000000000000;0.545300000000000;2.02600000000000;0.000100000000000000;0.800000000000000;0.344200000000000]
  
  
  
  %% trial 1
  
ymin = - model.G.Ly / 2 ;
ymax = + model.G.Ly / 2 ;

ypos = linspace( ymin , ymax , 89) ;
cummlativePowerAbsor = zeros(size(ypos)) ;
for j = 1 : length(ypos)

    model.MC.beam.yFocus = ypos(j) ; % [cm] y position of focus
    model.MC.LC.y        = ypos(j) ; % [cm] y position

    % Execution, do not modify the next line:
    model = runMonteCarlo(model);

    % plot(model,'MC');


    reflectance_vessel
    cummlativePowerAbsor(j)  = sum(cum_ref) ;

end

figure,
    plot(ypos , cummlativePowerAbsor)