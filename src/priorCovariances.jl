
#Converts 1-\sigma errors and a correlation matrix to a covariance matrix
function createPriorCovarianceMatrix(σError, correlationMatrix )
    correlationMatrix .* (σError * σError')
end

# Create a correlation matrix with correlation length scale (cls, in p coordinates), pressure grid p, and pressure at the UTLS (as it decorrelates Trop and Strat) 
function createCorrelationMatrix(p,p_utls,cls)
    correlationMatrix = zeros(length(p), length(p));
    for i=1:length(p)
        for j=1:length(p)
            if i==j
                correlationMatrix[i,j] = 1
            else
                if (p[i] < p_utls && p[j] > p_utls) || (p[i] > p_utls && p[j] < p_utls)
                    correlationMatrix[i,j] = 0.0;
                else
                    correlationMatrix[i,j] = exp(-abs(p[i]-p[j])/cls)
                end
            end
        end
    end
    return correlationMatrix
end

# This needs to be refined with different options, for now we assume higher error in the BL, 
# decrease towards the UTLS and known in the strat (can later add co-variances in the N2O and CH4 strat). 
function createErrorVector(p,p_utls,p_bl, rel_Error, vmr; bl_error=5, utls_error=0.0001)
    n_layers = length(p)
    σ = zeros(n_layers)
    for i=1:n_layers
        if p[i] > p_bl
            σ[i] = rel_Error * bl_error * vmr[i]
        elseif p[i] < p_utls
            σ[i] = rel_Error * utls_error * vmr[i]
        else
            σ[i] = rel_Error * vmr[i]
        end
    end
    return σ
end