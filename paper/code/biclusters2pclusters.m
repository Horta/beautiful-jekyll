function pclusters = biclusters2pclusters(biclusters, nr, nc) %#ok<INUSD>

%nr: number of rows
%nc: number of columns

k = length(biclusters);
pclusters = cell(1,k);

for i = 1:k
    
    rows = biclusters(i).rows;
    cols = biclusters(i).cols;
    
    pclusters{i} = [];
    for s = 1:length(cols)
        
        pclusters{i} = [pclusters{i} rows + nr * (cols(s)-1)];
    end
end

end