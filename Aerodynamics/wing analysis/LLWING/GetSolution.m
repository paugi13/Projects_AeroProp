function [cl_local, c4nods, force_coeff, chord] = GetSolution(N, ALPHA, FlapCorr, ...
    YF_pos, CF_ratio, DE_flap, A0p, CM0p, CDP, AR, TR, DE25, ETIP)

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = ...
    geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using plane Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,...
    chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

end

