#compute (log10_kappa, max_gamma_eager, max_gamma_dcr) for symmetric 6j
function compute_amplification(j, k::Int)
    h = k + 2
    
    # Symmetric 6j-symbol bounds
    z_min, z_max = 3j, 4j
    
    # Precompute  BigFloat q-integers
    qint_big = zeros(BigFloat, z_max + 1)
    log_qint = zeros(Float64, z_max + 1) # log[n]
    
    big_pi = big(π)
    for n in 1:z_max+1
        val = sin(n * big_pi / h) / sin(big_pi / h)
        qint_big[n] = abs(val)
        log_qint[n] = Float64(log(qint_big[n]))
    end

    # Möbius inversion for log|Phi_d(q^2)|
    log_Phi = zeros(Float64, z_max + 1)
    for n in 2:z_max+1
        val = log_qint[n]
        for d in 2:n-1
            if n % d == 0
                val -= log_Phi[d]
            end
        end
        log_Phi[n] = val
    end

    #track summation
    max_gamma_eager = 0.0
    max_gamma_dcr   = 0.0
    sum_S = big(0.0)
    sum_abs_T = big(0.0)
    
    #loop
    for z in z_min:z_max
        
        # compute kappa (κ)
        T_num = prod(qint_big[1:z+1])
        T_den = prod(qint_big[1:(z-3j)])^4 * prod(qint_big[1:(4j-z)])^3
        T_z = T_num / T_den
        
        sum_abs_T += T_z
        sum_S += (iseven(z) ? T_z : -T_z)
        
        # eager γ (Log-Sum swell)
        log_Nz_eager = sum(view(log_qint, 1:z+1))
        log_Dz_eager = 4 * sum(view(log_qint, 1:z-3j)) + 3 * sum(view(log_qint, 1:4j-z))
        
        gamma_eager = log_Nz_eager + log_Dz_eager
        (gamma_eager > max_gamma_eager) && (max_gamma_eager = gamma_eager)
        
        # --- C. DCR Gamma (Symbolic minimal swell) ---
        log_Nz_dcr, log_Dz_dcr = 0.0, 0.0
        
        for d in 2:z+1
            # Exact integer exponent sieve
            e_d = floor(Int, (z+1)/d) - 4 * floor(Int, (z-3j)/d) - 3 * floor(Int, (4j-z)/d)
            
            if e_d > 0
                log_Nz_dcr += e_d * log_Phi[d]
            elseif e_d < 0
                log_Dz_dcr -= e_d * log_Phi[d] # -e_d handles the absolute value
            end
        end
        
        gamma_dcr = log_Nz_dcr + log_Dz_dcr
        (gamma_dcr > max_gamma_dcr) && (max_gamma_dcr = gamma_dcr)
    end
    
    # Convert κ to log10
    log10_kappa = Float64(log10(sum_abs_T / abs(sum_S)))
    
    #convert gammas to log10 
    max_gamma_eager /= log(10)
    max_gamma_dcr /= log(10)
    
    return (log10_kappa, max_gamma_eager, max_gamma_dcr)
end

# j= 10
# compute_amplification(j,4j)
# (1.2748894577139092, 61.091413429608444, 18.996838137590856)

# j= 500
# compute_amplification(j,4j)
# (72.70031873685633, 9518.545434913727, 1079.7909965294405)