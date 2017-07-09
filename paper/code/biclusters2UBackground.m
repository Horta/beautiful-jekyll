function U = biclusters2UBackground(clusters, nr, nc, retsparse)

%nr: number of rows
%nc: number of columns

k = length(clusters);

noise = ones(1,nr*nc);

for i = 1:k
    
    rows = clusters(i).rows;
    cols = clusters(i).cols;
    
    for s = 1:length(cols)
        
        noise(rows + nr * (cols(s)-1)) = 0;
    end
end

nnoise = sum(noise);

U = zeros(k+nnoise,nr*nc);

for i = 1:k
    
    rows = clusters(i).rows;
    cols = clusters(i).cols;
    
    for s = 1:length(cols)
        
        U(i,rows + nr * (cols(s)-1)) = 1;
    end
end

noise = find(noise==1);
t = k+1:k+nnoise;
U(t + ((noise-1) * (k+nnoise))) = 1;


if nargin == 4 && retsparse
    
    U = sparse(U);
end

end