 function [q,r] = cGram_Schmidt(A)
% Computes the QR factorization of $A$ via classical Gram Schmid 

[~,m] = size(A);
q = A;
for j = 1:m
    for i = 1:j-1
        r(i,j) = q(:,j)'*q(:,i);
    end
    for i = 1:j-1
        q(:,j) = q(:,j)-r(i,j)*q(:,i);
    end
    t = norm(q(:,j),2);
    q(:,j) = q(:,j)/t;
    r(j,j) = t;
end