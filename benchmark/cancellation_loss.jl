
# Evaluate the catastrophic cancellation (Δ_loss) for the symmetric 6j-symbols
# compute digits loss

function cancellation_loss(j, k::Int)
    h = k + 2
    z_min, z_max = 3j, 4j
    
    # Precompute qints
    qint_big = zeros(BigFloat, z_max + 1)
    big_pi = big(π)
    for n in 1:z_max+1
        val = sin(n * big_pi / h) / sin(big_pi / h)
        qint_big[n] = abs(val)
    end

    # Compute prefactor 
    j_fact = prod(view(qint_big, 1:j))
    tri_fact = prod(view(qint_big, 1:3j+1))
    prefactor = ( (j_fact^3) / tri_fact )^2

    max_term_mag = big(0.0)
    sum_S = big(0.0)
    
    # loop
    for z in z_min:z_max
        T_num = prod(view(qint_big, 1:z+1))
        T_den = prod(view(qint_big, 1:z-3j))^4 * prod(view(qint_big, 1:4j-z))^3
        T_z = T_num / T_den
        
        if T_z > max_term_mag
            max_term_mag = T_z
        end
        
        sum_S += (iseven(z) ? T_z : -T_z)
    end
    
    # get true magnitude
    true_max_term = prefactor * max_term_mag
    final_mag = prefactor * abs(sum_S)
    
    # Δ_loss (Digits of precision destroyed)
    # even though the prefactor cancels out in the division, we apply it 
    # above so the returned magnitudes match the real physical values.
    delta_loss = Float64(log10(max_term_mag / abs(sum_S)))
    
    return (Float64(true_max_term), Float64(final_mag), delta_loss)
end



# j=50
# cancellation_loss(j,4j)
# (6031.4308321784365, 0.002188311033037061, 6.440311301754783)


# j=100
# cancellation_loss(j,4j)
# (2.9642274685145355e10, 0.0007795244429423376, 13.580081789848354)