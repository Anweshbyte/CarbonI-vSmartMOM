using SpecialPolynomials#, SparseArrays
export forward_model_x_ 

function forward_model_x(𝐱::AbstractArray{FT} ;sun = solarIrr, instrument=lociBox, sza=sza, vza = 0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	#@show dims
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show poly
    # Air Mass Factor
    AMF = 1/cosd(sza);# + 1/cosd(vza)
    
    # Total sum of τ
    ∑τ = zeros(FT,size(σ_matrix,1))
	#@show size(vmrs,2)
    for i=1:size(vmrs,2)
        #@show i, vmrs[end,i]
        #@show vmrs[:,i]
         ∑τ[:] += sum(σ_matrix[:,:,i] .* (vmrs[:,i] .* profile.vcd_dry)', dims=2)
         #@show i, maximum(∑τ), maximum(σ_matrix[:,:,i])
    end
    # Transmission without Tsolar
    T = sun .* reverse(exp.(-AMF * ∑τ))
	#@show ∑τ
    @time T_conv = CarbonI.conv_spectra(instrument, wl, T)
    L = T_conv;
    # x-axis for polynomial [-1,1], enables legendre later:
    x_poly = CarbonI.rescale_x(instrument.ν_out)
   return L .* poly.(x_poly) 
end

function forward_model_sat_x(𝐱::AbstractArray{FT} ; sun = solarIrr, instrument=wl_ci, convMatrix=cM, sza=sza, vza = 0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	#@show dims
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show poly
    # Air Mass Factor
    AMF = 1/cosd(sza) + 1/cosd(vza);
    
    # Total sum of τ
    ∑τ = zeros(FT,size(σ_matrix,1))
	#@show size(vmrs,2)
    for i=1:size(vmrs,2)
        for j=1:dims[2]
            @views ∑τ[:] .+= σ_matrix[:,j,i] .* (vmrs[j,i] .* profile.vcd_dry[j])
        end
    end
    # Transmission without Tsolar
    T = sun .* reverse(exp.(-AMF * ∑τ))
    
    T_conv = cM * T
    #L = T_conv;
    # x-axis for polynomial [-1,1], enables legendre later:
    x_poly = CarbonI.rescale_x(instrument)
    #x_poly = CarbonI.rescale_x(wl)
   return T_conv .* poly.(x_poly)
end

function forward_model_sat_x2(𝐱::AbstractArray{FT} ; sun = solarIrr::Vector{Float64}, indHR=indHR,indHR2=indHR2, indLR = indLR::UnitRange{Int64}, instrument=wl_ci::Vector{Float64}, convMatrix=cM::Matrix{Float64}, sza=sza::Float64, vza = 0.0, vcd_dry=profile.vcd_dry::Vector{Float64},σ_matrix=σ_matrix::Array{Float64,3})  where{FT}
    dims = size(σ_matrix)
	#@show dims
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show poly
    # Air Mass Factor
    AMF = 1/cosd(sza) + 1/cosd(vza);
    
    # Total sum of τ
    ∑τ = zeros(FT,length(indHR))
	#@show size(vmrs,2)
    for i=1:size(vmrs,2)
        for j=1:dims[2]
            @views ∑τ[:] .+= σ_matrix[indHR,j,i] .* (vmrs[j,i] .* vcd_dry[j])
        end
    end
    # Transmission without Tsolar
    @views T = reverse(exp.(-AMF * ∑τ))
    
    @views T_conv = cM[indLR,indHR2] * T
    #L = T_conv;
    # x-axis for polynomial [-1,1], enables legendre later:
    @views x_poly = CarbonI.rescale_x(instrument[indLR])
    #x_poly = CarbonI.rescale_x(wl)
   return T_conv .* poly.(x_poly)
end

function lin_forward_model_sat_x(𝐱::AbstractArray{FT} ;∑τ=∑τ, sun = solarIrr, instrument=lociBox, convMatrix=cM, sza=sza, vza = 0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	#@show dims
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show poly
    # Air Mass Factor
    AMF = 1/cosd(sza) + 1/cosd(vza);
    
    # Total sum of τ
    ∑τ = zeros(FT,size(σ_matrix,1))
    # Jacobi Matrix K
    K  = zeros(FT,length(lociBox.ν_out),length(𝐱));

    for i=1:size(vmrs,2)
        for j=1:dims[2]
            @views ∑τ[:] .+= σ_matrix[:,j,i] .* (vmrs[j,i] .* profile.vcd_dry[j])
        end
    end
    # Transmission without Tsolar
    T = sun .* reverse(exp.(-AMF * ∑τ))
    @time x_poly = CarbonI.rescale_x(1:length(sun))
    # Convolve and resample:
    L  = T .* poly.(x_poly)
    # x-axis for polynomial [-1,1], enables legendre later:
    
    # Multiply with Polynomial
    L_conv  = cM * L

    # Jacobian wrt to trace gases:
    for i=1:size(vmrs,2)
        for j=1:dims[2]
            d = AMF .* L .* -reverse(σ_matrix[:,j,i]) .* profile.vcd_dry[j]
            #how size(d)
            @views K[:,(i-1)*dims[2] + j] .= cM * d
        end
    end
    
    #x_poly = CarbonI.rescale_x(wl)
   return L_conv, K
end

function forward_model_x_(𝐱::AbstractArray{FT} ;sun = solarIrr,reflectance=refl, instrument=lociBox, sza=sza, vza=0.0, profile=profile,σ_matrix=σ_matrix, wl=wl) where {FT}
    dims = size(σ_matrix)
	# @show dims
    #xx
    vmrs = reshape(𝐱[1:(dims[2]*dims[3])],(dims[2],dims[3]) )
    poly = Legendre(𝐱[dims[2]*dims[3]+1:end])
    #@show size(vmrs)
    # Air Mass Factor
    AMF = 1/cosd(sza) + 1/cosd(vza)
    #@show sza
    # Total sum of τ
    ∑τ = zeros(FT,size(σ_matrix,1))
	#@show size(vmrs,2)
    for i=1:size(vmrs,2)
        #@show i, vmrs[end,i]
         ∑τ[:] += sum(σ_matrix[:,:,i] .* (vmrs[:,i] .* profile.vcd_dry)', dims=2)
    end
    # Transmission without Tsolar
    T = sun .* reflectance .* reverse(exp.(-AMF * ∑τ))
	#@show T
    T_conv = CarbonI.conv_spectra(instrument, wl, T)
    L = cosd(sza)*T_conv/π;
    # x-axis for polynomial [-1,1], enables legendre later:
    x_poly = CarbonI.rescale_x(instrument.ν_out)
    #@show poly.(x_poly)
   return L .* poly.(x_poly)
end