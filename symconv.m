function [c] = symconv(f,g,k)
   N = length(f)-1;
   shift = N/2+1;
   c = sym(zeros(N+1,1));
   for k=-N/2:N/2
     for p=-N/2:N/2
       for q=-N/2:N/2
         if ( (p + q) == k)
           c(k + shift,1) = c(k+shift,1) + f(p + shift)*g(q + shift);
         end
       end
     end
   end
end 
