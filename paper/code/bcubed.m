function r = bcubed(U, V)

% A. Bagga and B. Baldwin, "Entity-based cross-document coreferencing
% using the vector space model,? in Proceedings of the 17th international
% conference on Computational linguistics - Volume 1, ser. COLING ?98.
% Stroudsburg, PA, USA: Association for Computational Linguistics, 1998,
% pp. 79?85.

Mu = U'*U;
Mv = V'*V;

prec = bcubedPrecision(Mu, Mv);
recall = bcubedRecall(Mu, Mv);

r = ( 0.5*(1/prec + 1/recall) )^(-1);

end

function r = bcubedPrecision(Mu, Mv)

n = size(Mu,1);

MM = Mu .* Mv;

r = 0;
for i = 1:n
    
    indices = Mu(i,:) > 0;
    
    prec = MM(i,indices) ./ Mu(i,indices);
    
    r = r + mean(prec);
end

r = r/n;

end

function r = bcubedRecall(Mu, Mv)

r = bcubedPrecision(Mv, Mu);

end