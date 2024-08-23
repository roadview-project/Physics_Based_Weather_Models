function C_back = C_BACK(x,m,N_Max)

    sum = 0;
    for n=1:N_Max
        psi_x = x*sqrt(pi/(2*x))*besselj(n+(1/2),x);
        psi_minusOne_x = x*sqrt(pi/(2*x))*besselj((n-1)+(1/2),x);
        psi_mx = m*x*sqrt(pi/(2*m*x))*besselj(n+(1/2),m*x);
        psi_minusOne_mx = m*x*sqrt(pi/(2*m*x))*besselj((n-1)+(1/2),m*x);
    
        rn_mx = psi_minusOne_mx/psi_mx;
    
        epsilon_x_Real = psi_x;
        epsilon_x_Imag = x*sqrt(pi/(2*x))*bessely(n+(1/2),x);
        epsilon_x = complex(epsilon_x_Real,epsilon_x_Imag);
        epsilon_minusOne_x_Real = psi_minusOne_x;
        epsilon_minusOne_x_Imag = x*sqrt(pi/(2*x))*bessely((n-1)+(1/2),x);
        epsilon_minusOne_x = complex(epsilon_minusOne_x_Real,epsilon_minusOne_x_Imag);
    
        numerator_an = ( ( (rn_mx/m)+(n*(1-(1/(m^2)))/x) ) * psi_x) - psi_minusOne_x;
        denominator_an = ( ( (rn_mx/m)+(n*(1-(1/(m^2)))/x) ) * epsilon_x)-epsilon_minusOne_x;
        an = numerator_an/denominator_an;
    
        numerator_bn = (rn_mx*m*psi_x)-psi_minusOne_x ;
        denominator_bn = (rn_mx*m*epsilon_x)-epsilon_minusOne_x;
        bn = numerator_bn/denominator_bn;
    
        sum = sum + ( (2*n + 1)*((-1)^n)*(an-bn) );
    end
    C_back = (real(sum)^2)+(imag(sum)^2);

end