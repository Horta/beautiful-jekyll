function [biclusters,scores] = bcca(X, theta, r)
% 2009, Bi-correlation clustering algorithm for determining
% a set of co-regulated genes
% theta: Pearson correlation threshold
% r: minimum number of columns

[n p] = size(X);

biclusters = [];

for i = 1:n-1
    for j = i+1:n
        
        I = [i j];
        J = 1:p;
        
        while corr(X(I(1),J)', X(I(2),J)') < theta && length(J) >= r
            
            imaior = find_max_increase(X(I(1),J), X(I(2),J));
            J = J([1:imaior-1 imaior+1:length(J)]);
        end
        
        if corr(X(I(1),J)', X(I(2),J)') >= theta && length(J) >= r
            
            genes = [1:i-1 i+1:j-1 j+1:n];
            for ii = 1:length(genes)
                
                yes = iscompatible(X(genes(ii),J), X, I, J, theta);
                if yes
                    I = [I genes(ii)];
                end
            end
            exist = check_existence(biclusters, I, J);
            if ~exist
                ii = length(biclusters)+1;
                biclusters(ii).rows = I;
                biclusters(ii).cols = J;
            end
        end
        
        
    end
end

scores = evaluate_biclusters(biclusters, X);

[scores indices] = sort(scores,'descend');

biclusters = biclusters(indices);

end

function scores = evaluate_biclusters(biclusters, X)

k = length(biclusters);
scores = zeros(1,k);

for i = 1:k
    
    scores(i) = sum(sum(triu(corr(X(biclusters(i).rows,biclusters(i).cols)'),1))) / nchoosek(length(biclusters(i).rows),2);
end

end

function exist = check_existence(biclusters, I, J)

exist = false;
I = sort(I);
J = sort(J);
for i = 1:length(biclusters)
    if length(biclusters(i).rows) == length(I) && all(biclusters(i).rows==I) &&...
            length(biclusters(i).cols) == length(J) && all(biclusters(i).cols==J)
        
        exist = true;
        break
    end
end
end

function yes = iscompatible(g, X, I, J, theta)

yes = true;
for i = 1:length(I)
    
    if corr(g', X(I(i),J)') < theta
        yes = false;
        break
    end
end

end


function imaior = find_max_increase(g1, g2)

p = length(g1);

maior = -inf;
imaior = -1;
for i = 1:p
    
    t = corr(g1([1:i-1 i+1:p])', g2([1:i-1 i+1:p])');
    if t > maior
        maior = t;
        imaior = i;
    end
end

assert(imaior ~= -1, 'Could not calculate the correlation between g1 and g2.');

end