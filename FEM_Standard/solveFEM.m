function U = solveFEM(A,F)

U = A\F;

U = [0;U;0];

end
