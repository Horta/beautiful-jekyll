function N = anne_num_cluster_belonging_to(biclusters, n, p)

N = zeros(n,p);

for i = 1:length(biclusters)
    
    cols = biclusters(i).cols;
    rows = biclusters(i).rows;
    
    for j = 1:length(cols)
        
        N(rows,cols(j)) = N(rows,cols(j)) + 1;
    end
end

end