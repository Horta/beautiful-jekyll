function clusters = U2clusters( U )

k = size(U,1);

clusters = cell(1,k);

for i = 1:k
 
    clusters{i} = find( U(i,:) == 1 );
end

end
