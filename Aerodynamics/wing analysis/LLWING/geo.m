function [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr)

S = 1.0/AR ;   % wing surface

cr = 2.0*S/(1.0+TR) ;  % root chord

h = 1.0 / N  ; % strip span

c4nods = zeros(3,N) ; % quarter chord nodes (panels' midpoints)
chord = zeros(N,1) ; % panel chord (at midpoint)
c75nods = zeros(3,N) ; % panel control points 
s_pan = zeros(N,1)  ; % panel surface area

c4nods(2,1) = -0.5 + h/2. ; 

c4nods(3,:) = 0.0 ;  % (just in case of dihedral, by the time being only planar wing)

for i = 2:N % Y position (bounded vortex's midpoint)
    c4nods(2,i) = c4nods(2,i-1) + h ; 
end

tande25 = tan(DE25*pi/180.) ;  

for i = 1:N 
    
    c4nods(1,i) = 0.25*cr + abs(c4nods(2,i))*tande25 ;  % x position (bounded vortex's midpoint)
    
    chord(i) = cr*(1.0-(1.-TR)*2.0*abs(c4nods(2,i))) ;   % panel's chord (at midpoint)
    
    c75nods(:,i) = c4nods(:,i) + [chord(i)*0.5,0.,0.].' ;  % control points
    
    s_pan(i) = chord(i)*h ;  % panel surface area
    
end

S = sum(s_pan) ;  % wing surface area (=1/AR)

% Airfoil flap effects (from thin airfoil theory)

xh = 1.0 - CF_ratio  ; th = acos(1-2*xh)  ;   % hinge position

eta = DE_flap*FlapCorr ; % corrected flap deflection (deg)

delta_a0_flap = -(1-th/pi + sin(th)/pi)*eta   ; % delta_alfa0 (in deg) 
delta_Cm0_flap = -0.5*sin(th)*(1-cos(th))*eta*pi/180 ; % delta_CM0

% Section values

normals = zeros(3,N) ; eps = zeros(N,1) ; % by the time being only planar wing
Cm0_y = zeros(N,1) ; a0_y = zeros(N,1) ;

cm0_root = CM0p(1) ; cm0_tip = CM0p(2) ; % sections' free moments
a0_root = A0p(1) ; a0_tip = A0p(2) ;  % sections' angles of zero lift

for i = 1:N
    
    ypos = abs(c4nods(2,i))*2. ;  % normalized span position
    
    eps(i) = ETIP*ypos ;  % geometric twist
    
    a0_y(i) = a0_root + (a0_tip-a0_root)*ypos ;  % section alpha zero-lift
    Cm0_y(i) = cm0_root + (cm0_tip-cm0_root)*ypos ;  % section free moment
    
    if ( ypos >= YF_pos(1) && ypos <= YF_pos(2) )  % flap incremental values
        a0_y(i) = a0_y(i) + delta_a0_flap ;
        Cm0_y(i) = Cm0_y(i) + delta_Cm0_flap ;
    end
    
    zzl_angle = -a0_y(i) + eps(i) ;  % section's total zero-lift line angle
    
    normals(1,i) = sind(zzl_angle);
    normals(2,i) = 0.0 ;
    normals(3,i) = cosd(zzl_angle);
    
end

mac = 2.0/3.0*cr*(1.+TR+TR^2)/(1.+TR) ;  % mean aerodynamic chord

end  % function