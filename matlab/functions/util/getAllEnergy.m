% get single and pairwise potential for all
function U = getAllEnergy(unary, pairwise, labels, neighbors)
N = size(pairwise, 1);
U = 0;
for i = 1:N
    neigh = neighbors{i};
    notSame = find(labels(neigh) ~= labels(i));
    U = U + unary(labels(i)+1, i) + sum(pairwise(i, neigh(notSame)));
end

function [fobj] = getAllEnergy2(unary, pairwise, labels, neighbors)
N = size(pairwise, 1);
fobj=0;
for i=1:N,
    fobj = fobj + unary(labels(i)+1, i);
    neigh = neighbors{i};
    v1=labels(neigh);
    for j=1:length(neigh),
        if(v1(j)~=labels(i))
            fobj = fobj + pairwise(i, neigh(j));
        end
    end
end


