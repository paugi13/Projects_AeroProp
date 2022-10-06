function [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals)

ncases = size(ALPHA,2) ; RHS = zeros(N,1) ; GAMMA = zeros(N,ncases) ;

Ui = zeros(3,ncases) ;

for icase = 1:ncases   % get circulations for different angles of attack
    
    aoa = ALPHA(icase) ;
    
    Ui(:,icase) = [ cosd(aoa),0.0,sind(aoa) ].' ;  % freestream velocity
    
    for i = 1:N    % right-hand side
        RHS(i) = -dot(Ui(:,icase),normals(:,i)) ;
    end
    
    GAMMA(:,icase) = inv_A*RHS(:) ;
    
end

end % function