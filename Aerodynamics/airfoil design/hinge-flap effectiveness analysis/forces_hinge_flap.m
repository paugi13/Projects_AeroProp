function [cl_hinge, cm_hinge, force_hinge, moment_hinge] = ...
    forces_hinge_flap(point, pan, f, p, t, alpha, xh, def)

%% FUNCTION TO CALCULATE HINGE CONTRIBUTIONS (COEFFICIENT VARIATIONS DUE TO FLAPS).
%Angle of attack.
u_inf = 100;                                      %Freestream.
dens = 1.225;                                       %Density.
c = 1;
pos = zeros(point, 2);

for i = 1:point
    pos(i, 1) = 1/2*(1-(cos((i-1)/(point-1)*pi)));
end

for i = 1:point
    if pos(i, 1)<=p
        pos(i,2) = f/(p^2)*(2*p*pos(i, 1)-pos(i, 1)^2);
        
    else
        pos(i,2) = f/(1-p)^2*(1-2*p+2*p*pos(i, 1)-pos(i, 1)^2);
    end
end

if p==0 && f==0
    pos(:,2)=0;
end

% if def ~= 0  
    
%% Find xh approx
d_min = 1000;
for i = 1:point
    if abs(pos(i,1)-xh)<d_min
        d_min = abs(pos(i,1)-xh);
        auxe = i;
    end
end


%% Rotate points behind xh
coord_rel = zeros(2,1);
for i = 1:point
    if pos(i, 1)>pos(auxe, 1)
        coord_rel(:, 1) = pos(i, :) - pos(auxe, :);
        rot_mat = [cosd(-def) -sind(-def); sind(-def) cosd(-def)];
        new_coord_rel = rot_mat*coord_rel;
        pos(i, 1) = pos(auxe, 1) + new_coord_rel(1);
        pos(i, 2) = pos(auxe, 2) + new_coord_rel(2);
    end
end
    
% end

%% Normal vector calculus
%Positive vector of the panel is also calculated
n = [0,0,1];
vec_tan = zeros(pan, 2);
vec_norm = zeros(pan, 2);
for i = 1:pan
    vec_tan(i, :) = pos(i+1,:) - pos(i,:);
    vec_incr = [vec_tan(i, 1), vec_tan(i, 2), 0];
    v_normal = cross(n, vec_incr);
    vec_norm(i, 1) = v_normal(1)/norm(v_normal);
    vec_norm(i, 2) = v_normal(2)/norm(v_normal);
end

%% Vortex and control point position
%One per panel 
%vec_tan is not unit vector
cp = zeros(pan, 2);
vortex_point = zeros(pan, 2);
for i = 1:pan
    cp(i, :) = pos(i, :) + 0.75*vec_tan(i, :);
    vortex_point(i, :) = pos(i, :) + 0.25*vec_tan(i, :);
end

%% Vortex equation system
%pan x pan system
%Effects on control points from vortex. 

A = zeros(pan, pan);
rhs = zeros(pan, 1);
for i = 1:pan
    for j = 1:pan
        vec_dist = [0 0];
        vec_dist(1) =  cp(i, 1) - vortex_point(j, 1);
        vec_dist(2) =  cp(i, 2) - vortex_point(j, 2);
        r = norm(vec_dist);
        u = 1/(2*pi)*vec_dist(2)/(r^2);
        w = -1/(2*pi)*vec_dist(1)/(r^2);
        %Every row is assigned to one panel. 
        %Every column is related to the effects of one specific vortex. (Ej:
        %column 1 is related to the effects of the vortex from the first
        %panel).
        A(i, j) = u*vec_norm(i, 1) + w*vec_norm(i,2);     
    end
    rhs(i) = -u_inf*(cosd(alpha)*vec_norm(i, 1) + sind(alpha)*vec_norm(i, 2));
end

circ_pan = A\rhs;   %System solver.

%% Cl and Cm calculus

size_1 = size(circ_pan);
sum = 0;

%Only the flaps lift. It starts adding values from the vortex at position
%aux.
for i = 1:size_1(1)
    if i>= auxe
    sum = sum + circ_pan(i);
    end
end

cl_hinge = 2*sum/(u_inf*c);
force_hinge = dens*sum*u_inf;

%Momentum around LE (x_ref=0)
%-2/u*c^2 sum(circ*)*cos al

sum = 0;
x_ref = pos(auxe,1);

for i = 1:pan   
    if i>= auxe
     incr_mom = circ_pan(i)*(vortex_point(i, 1)-x_ref);
     sum = sum + incr_mom;
    end
end
cm_hinge = -2/(u_inf*c^2)*sum*cosd(alpha);
moment_hinge = cm_hinge*dens*0.5*(u_inf^2)*(c^2);

end