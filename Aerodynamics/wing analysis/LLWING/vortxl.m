
% ·································
% Velocity induced by a vortex line
% ·································

function [u] = vortxl (X1,X2,XP,gamma)

inv_4pi = 0.25/pi;

r0 = X2-X1 ;
r1 = XP-X1 ;  % defined as -(X-XP) to have beta in the right side
r2 = XP-X2 ;

norm_r1 = norm(r1) ;
norm_r2 = norm(r2) ;

r1xr2 = cross(r1,r2) ;
norm_r1xr2 = norm(r1xr2) ;

tol = 1.0e-6 ; u = zeros(3,1) ;

if ( norm_r1>tol && norm_r2>tol && norm_r1xr2>tol )
    
    inv_r1xr2 = 1.0 / norm_r1xr2 ;
    inv_r1 = 1.0 / norm_r1 ;
    inv_r2 = 1.0 / norm_r2 ;
    
    a = r0 * inv_r1xr2 ;
    b = r1*inv_r1 - r2*inv_r2 ;
    c = dot(a,b) ;
    
    u = inv_4pi*gamma*c*r1xr2*inv_r1xr2 ;
    
else
    u = 0.0 ;
end

end % function