function [cl, cm_le, pos] = cl_cm0_2408_flap_function(point, pan, f, p, t, alpha, xh, def)
                                 %Angle of attack.
u_inf = 1;                                      %Freestream.
dens = 1;                                       %Density.
c = 1;

i=1;
pos = zeros(point, 2);
while i<=point
    pos(i, 1) = 1/2*(1-(cos((i-1)/(point-1)*pi)));     %Cosine dist.
%     pos(i,1) = c*(i-1)/(point-1);                       %Linear dist.
    i=i+1;
end

i=1;

while i<=point
   
    if pos(i, 1)<=p
        pos(i,2) = f/(p^2)*(2*p*pos(i, 1)-pos(i, 1)^2);
        
    else
        pos(i,2) = f/(1-p)^2*(1-2*p+2*p*pos(i, 1)-pos(i, 1)^2);
    end
    i = i+1;
end

%% Define xh and coordinates

% if xh <= p 
%     zh = f/(p^2)*(2*p*xh-xh^2);
% else
%     zh = f/(1-p)^2*(1-2*p+2*p*xh-xh^2);
% end

%% Find xh approx

i = 1;
d_min = 1000;
while i<=point
    if abs(pos(i,1)-xh)<d_min
        d_min = abs(pos(i,1)-xh);
        aux = i;
    end
    i = i+1;
end


%% Rotate points behind xh

i = 1;
coord_rel = zeros(2,1);
while i<=point
    if pos(i, 1)>pos(aux, 1)
        coord_rel(:, 1) = pos(i, :) - pos(aux, :);
        rot_mat = [cosd(-def) -sind(-def); sind(-def) cosd(-def)];
        new_coord_rel = rot_mat*coord_rel;
        pos(i, 1) = pos(aux, 1) + new_coord_rel(1);
        pos(i, 2) = pos(aux, 2) + new_coord_rel(2);
    end
    i=i+1;
end

%% Define thickness
%Here we could define the thickness (optional)

%% Normal vector calculus
%Positive vector of the panel is also calculated
i = 1;
n = [0,0,1];
vec_tan = zeros(pan, 2);
vec_norm = zeros(pan, 2);
while i<=pan
    
    vec_tan(i, 1) = pos(i+1,1) - pos(i,1);
    vec_tan(i, 2) = pos(i+1,2) - pos(i,2);
    vec_incr = [vec_tan(i, 1), vec_tan(i, 2), 0];
    v_normal = cross(n, vec_incr);
    vec_norm(i, 1) = v_normal(1)/norm(v_normal);
    vec_norm(i, 2) = v_normal(2)/norm(v_normal);
    
     i = i+1;
end

%Check normal vectors.
% figure
% quiver(pos()vec_norm(3,1), vec_norm(3,2));

%% Vortex and control point position
%One per panel 
%vec_tan is not unit vector

i = 1;
cp = zeros(pan, 2);
vortex_point = zeros(pan, 2);
while i<=pan
    cp(i, 1) = pos(i, 1) + 0.75*vec_tan(i, 1);
    cp(i, 2) = pos(i, 2) + 0.75*vec_tan(i, 2);
    vortex_point(i, 1) = pos(i, 1) + 0.25*vec_tan(i, 1);
    vortex_point(i, 2) = pos(i, 2) + 0.25*vec_tan(i, 2);
    i = i+1;
end

%Check points
% figure
% plot(pos(:,1), pos(:,2));
% axis equal
% hold on 
% scatter(cp(:,1), cp(:,2));
% scatter(vortex_point(:,1), vortex_point(:,2))
% hold off


%% Vortex equation system
%pan x pan system
%Effects on control points from vortex. 

i=1;
j=1;

A = zeros(pan, pan);
rhs = zeros(pan, 1);

while i<=pan
    while j<=pan
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
        j=j+1;
    end
    rhs(i) = -u_inf*(cosd(alpha)*vec_norm(i, 1) + sind(alpha)*vec_norm(i, 2));
    j=1;
    i=i+1;
end

circ_pan = A\rhs;   %System solver.

%% Cl and Cm calculus

i = 1;
size_1 = size(circ_pan);
sum = 0;

while i<=size_1(1)
    sum = sum + circ_pan(i);
    i = i+1;
end

cl = 2*sum/(u_inf*c);

%Momentum around LE (x_ref=0)
%-2/u*c^2 sum(circ*)*cos al
i = 1;
sum = 0;
x_ref = 0;

while i<= pan
   incr_mom = circ_pan(i)*(vortex_point(i, 1)-x_ref);
   sum = sum + incr_mom;
   i = i+1;
end

cm_le = -2/(u_inf*c^2)*sum*cosd(alpha);


end
