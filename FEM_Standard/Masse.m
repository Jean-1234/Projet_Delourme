function M=Masse(x,h)
  n=length(x)-2;
  M=sparse(n+2,n+2);
  for k=1:n+1
    M(k:k+1,k:k+1)=M(k:k+1,k:k+1)+h*[1/3 1/6 ; 1/6 1/3] ;
  end
end
