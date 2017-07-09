function r = exbcubed(U, V)

% E. Amigo, J. Gonzalo, J. Artiles, and F. Verdejo, "A comparison of
% extrinsic clustering evaluation metrics based on formal constraints,"
% Inf. Retr., vol. 12, no. 4, pp. 461?486, Aug. 2009.

Mu = U'*U;
Mv = V'*V;

prec = ExBCubedPrecision(Mu, Mv);
recall = ExBCubedRecall(Mu, Mv);

r = ( 0.5*(1/prec + 1/recall) )^(-1);

end

function r = ExBCubedPrecision(Mu, Mv)

n = size(Mu,1);

MM = min(Mu, Mv);

r = 0;
for i = 1:n
    
    indices = Mu(i,:) > 0;
    
    prec = MM(i,indices) ./ Mu(i,indices);
    
    r = r + mean(prec);
end

r = r/n;

end

function r = ExBCubedRecall(Mu, Mv)

r = ExBCubedPrecision(Mv, Mu);

end