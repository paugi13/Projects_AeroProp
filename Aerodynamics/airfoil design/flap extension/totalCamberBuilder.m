function posGen = totalCamberBuilder(pos1,pos2)
% Function to build assemble a general node matrix. 

if size(pos1, 2) ~= size(pos2, 2)
    error('Input matrix do not have de same size.');
end

posGen = zeros(size(pos1,1)+size(pos2, 1), size(pos1, 2));

for i=1:size(pos1, 1)
    posGen(i, :) = pos1(i,:);
end

for i=1:size(pos2, 1)
    posGen(size(pos1, 1)+i, :) = pos2(i,:);
end
end

