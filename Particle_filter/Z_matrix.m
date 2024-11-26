function z=Z_matrix(w,signal)


for j=1:length(signal)
   
  z(j,1)= exp(1i*w*(j-1)); 
  z(j,2)= exp(2*1i*w*(j-1));
       
end

end