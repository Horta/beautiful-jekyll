function U = clusters2U(clusters)

n = max([clusters{:}]);
k = length(clusters);

U = zeros(k,n);

for i = 1:k
    
    U(i,clusters{i}) = 1;
end

end