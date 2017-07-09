function r = biclusteringError(biclustersa, biclustersb, n, p)

% the actual name is clustering error
% mas estou usando biclustering error porque eu especifiquei a implementacao
% para biclustering.

ka = length(biclustersa);
kb = length(biclustersb);

pa = biclusters2pclusters(biclustersa, n, p);
pb = biclusters2pclusters(biclustersb, n, p);

M = zeros(ka,kb);
for i = 1:ka
    for j = 1:kb
        
        M(i,j) = length( intersect(pa{i},pb{j}) );
    end
end

if ka < kb
    
    M = [M; zeros(kb-ka,kb)];
elseif ka > kb
    
    M = [M zeros(ka, ka-kb)];
end

[~,cost] = munkres(-M);

dmax = -cost;

Us = anne_union_size(biclustersa, biclustersb, n, p);

r = 1 - ((Us - dmax) / Us);
    
end
