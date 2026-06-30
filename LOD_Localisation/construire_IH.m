function IH = construire_IH(Xh, XH)

NH = length(XH) - 1;
Nh = length(Xh) - 1;

IH = sparse(NH+1, Nh+1);

tol = 1e-12;

[tf, loc] = ismembertol(XH, Xh, tol);

lignes = find(tf);
colonnes = loc(tf);

IH = sparse(lignes, colonnes, 1, NH+1, Nh+1);

end
