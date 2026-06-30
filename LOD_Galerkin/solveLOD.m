function U_lod = solveLOD(A_lod, F_lod)
% solveLOD : resout le systeme LOD et ajoute les CL homogenes.

U_lod = [0; A_lod \ F_lod; 0];

end
