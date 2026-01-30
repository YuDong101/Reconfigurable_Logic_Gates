function xm = find_xm_only(I, x_range, step)
    N = 3;
    a = zeros(1, N);
    xm = 0; 
    Um = 0;
    x = -x_range;
       
    while x <= x_range
        c=0.8;
        FI = 4*2*c * (1-I.^2) * I.^2 -c;
        U = 5*x^4 - 2*x^2 - FI*x;
       
        for i = 1:(N-1)
            a(i) = a(i+1);
        end
        a(N) = U;
        
        if x > ( - x_range) && a(2) > a(1) && a(2) > a(3)
                xm = x - step;
                Um = a(2);
        end
        x = x + step;
    end     
    
    if  abs(abs(xm) - 1.5) < 1e-6
          xm=0.0;
    end
% fprintf('xₘ = %.6f, Uₘ = %.6f\n', xm, Um);
end