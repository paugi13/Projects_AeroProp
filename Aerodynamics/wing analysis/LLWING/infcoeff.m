function [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
    
A = zeros(N,N) ; one = 1.0 ;

wake_len = 20. ; % assumes a 20b wake length

for i = 1:N
    
    xp = c75nods(:,i) ; np = normals(:,i) ; % control point
    
    for j = 1:N
        
        xb = c4nods(:,j) + [0.0, -0.5*h, 0.0].' ;  % horseshoe vortex (planar wake)
        xc = c4nods(:,j) + [0.0, 0.5*h, 0.0].' ;   
        xa = xb + [wake_len,0.,0.].';   
        xd = xc + [wake_len,0.,0.].';
        
        u1 = vortxl(xa,xb,xp,one) ;  % first trailing vortex
        u2 = vortxl(xb,xc,xp,one) ;  % bounded vortex
        u3 = vortxl(xc,xd,xp,one) ;  % second trailing vortex
        
        u_tot = u1+u2+u3 ; % induced velocity
        
        A(i,j) = dot(u_tot,np) ; % influence coefficient
        
    end
    
end

inv_A = inv(A) ;  % inverse of the coefficients matrix

end  % function