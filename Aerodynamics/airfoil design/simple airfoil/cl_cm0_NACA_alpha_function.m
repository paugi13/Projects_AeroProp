function [cl, cm_le] = cl_cm0_NACA_alpha_function(point, pan, alpha, pos, xh, def)
                                
u_inf = 1;                                                                           
c = 1;

if def ~= 0  
    

    d_min = 1000;
    for i = 1:point
        if abs(pos(i,1)-xh)<d_min
            d_min = abs(pos(i,1)-xh);
            aux = i;
        end
    end

    %% Rotate points behind xh

    coord_rel = zeros(2,1);
    for i = 1:point
        if pos(i, 1)>pos(aux, 1)
            coord_rel(1, 1) = pos(i, 1) - pos(aux, 1);
            coord_rel(2, 1) = pos(i, 2) - pos(aux, 2);
            rot_mat = [cosd(-def) -sind(-def); sind(-def) cosd(-def)];
            new_coord_rel = rot_mat*coord_rel;
            pos(i, 1) = pos(aux, 1) + new_coord_rel(1);
            pos(i, 2) = pos(aux, 2) + new_coord_rel(2);
        end
    end
    
end



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

%Check normal vectors.
% figure
% quiver(pos()vec_norm(3,1), vec_norm(3,2));

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
        vec_dist(1) = -vortex_point(j, 1) + cp(i, 1);
        vec_dist(2) = -vortex_point(j, 2) + cp(i, 2);
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

for i = 1:size_1(1)
    sum = sum + circ_pan(i);
end

cl = 2*sum/(u_inf*c);

%Momentum around LE (x_ref=0)
%-2/u*c^2 sum(circ*)*cos al
sum = 0;
x_ref = 0;

for i = 1:pan
   incr_mom = circ_pan(i)*(vortex_point(i, 1)-x_ref);
   sum = sum + incr_mom;
end

cm_le = -2/(u_inf*c^2)*sum*cosd(alpha);
end