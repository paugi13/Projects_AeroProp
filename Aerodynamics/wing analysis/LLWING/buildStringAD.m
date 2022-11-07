function numSt = buildStringAD(wantedAoA)

if mod(wantedAoA, 1) ~= 0
    if mod(wantedAoA*10, 1) ~= 0    % two decimals
        st1 = num2str(fix(wantedAoA));
        st2 = num2str(fix(mod(wantedAoA,1)*10));
        st3 = num2str(mod(mod(wantedAoA,1)*10, 1)*10);
        numSt = join([st1, st2, st3]);
    else % one decimal
        st1 = num2str(fix(wantedAoA));
        st2 = num2str(mod(wantedAoA, 1)*10);
        numSt = join([st1, st2]);
    end
else        % no remainder
    numSt = num2str(wantedAoA);
end

end

