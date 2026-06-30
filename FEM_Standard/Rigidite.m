function K=Rigidite(x,h)

  n=length(x)-2;
  K=sparse(n+2,n+2);

  for k=1:n+1
    K(k:k+1,k:k+1)=K(k:k+1,k:k+1)+ (1/h) *[1 -1 ; -1 1] ;
  end


end

