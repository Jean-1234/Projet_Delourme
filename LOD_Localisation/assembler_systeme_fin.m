function [Kh, Mh] = assembler_systeme_fin(Xh, a_eps)

    Nh = length(Xh) - 1;

    Kh = sparse(Nh+1, Nh+1);
    Mh = sparse(Nh+1, Nh+1);

    h = Xh(2:end) - Xh(1:end-1);
    x_mid = 0.5 * (Xh(2:end) + Xh(1:end-1));
    a_val = a_eps(x_mid);

    for j = 1:Nh

        idx = j:j+1;

        Kh(idx,idx) = Kh(idx,idx) + (a_val(j)/h(j)) * [1 -1; -1 1];

        Mh(idx,idx) = Mh(idx,idx) + (h(j)/6) * [2 1; 1 2];

    end

end
