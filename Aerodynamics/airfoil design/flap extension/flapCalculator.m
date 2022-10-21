function posFlap = flapCalculator(rotAngle, posFlap, pointFlap, xRel)
% Function that rotates the flap's camber line. 

coord_rel = zeros(2,1);
for i=1:pointFlap
    coord_rel(:, 1) = posFlap(i, :);
    rot_mat = [cosd(-rotAngle) -sind(-rotAngle); sind(-rotAngle)...
        cosd(-rotAngle)];
    new_coord_rel = rot_mat*coord_rel;
    posFlap(i, 1) = new_coord_rel(1) + xRel;
    posFlap(i, 2) = new_coord_rel(2);
end

end

