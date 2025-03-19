
using Unitful
using CarbonI
using Base

Base.@kwdef mutable struct InstrumentSpecs
    ### Constants; these are design characteristics
    # not expected to change
    ET = 44.0u"ms";
    Pitch = 18.0u"μm";
    Fnumber::Float64 = 2.2;
    coadd_rate::Int = 10;

    ### Dependent on build outcome; each has a 
    # required value and a current best estimate (CBE)
    FPA_quantum_efficiency::Float64 
    bench_efficiency::Float64
    readout_noise::Float64
    dark_current
    SSI
    FWHM::Float64
    lower_wavelength::Float64
    upper_wavelength::Float64 
    modelling_Δwl::Float64
    instrument_wl 
    modelling_wl
    convolution_matrix::Matrix{Float64}
    instrument_kernel::CarbonI.KernelInstrument
    pixel_size_global
    pixel_size_target

end

function build_instrument(;FPA_quantum_efficiency::Float64, 
                              bench_efficiency::Float64, 
                              readout_noise::Float64, 
                              dark_current,
                              SSI,
                              FWHM::Float64,
                              lower_wavelength::Float64,
                              upper_wavelength::Float64,
                              pixel_size_global,
                              pixel_size_target)

        modelling_Δwl = 0.005
        instrument_wl = collect(lower_wavelength:ustrip(SSI):upper_wavelength)
        if lower_wavelength - 5*ustrip(SSI) < 2030
            error("Instrument lower wavelength is below 2030 nm")
        end
        if upper_wavelength + 5*ustrip(SSI) > 2385
            error("Instrument upper wavelength is above 2385 nm")
        end
        # Stick to a common grid to make life easier for comparisons 
        modelling_wl = 2030:modelling_Δwl:2385
       
        cM, wl_ci, lociBox = CarbonI.create_carbonI_conv_matrix(modelling_wl,
	                                                            instrument_wl;
												                FWHM,
												                SSI,
												                verbose_return=true); 
       
        return InstrumentSpecs(FPA_quantum_efficiency=FPA_quantum_efficiency, 
                   bench_efficiency=bench_efficiency, 
                   readout_noise=readout_noise, 
                   dark_current=dark_current,
                   SSI=SSI,
                   FWHM=FWHM,
                   lower_wavelength=lower_wavelength,
                   upper_wavelength=upper_wavelength,
                   modelling_Δwl=modelling_Δwl,
                   instrument_wl=instrument_wl,
                   modelling_wl=modelling_wl,
                   convolution_matrix=cM,
                   instrument_kernel=lociBox,
                   pixel_size_global=pixel_size_global,
                   pixel_size_target=pixel_size_target)
end

function build_instrument(type::String)
    if type == "CBE"
        return cbe_instrument()
    elseif type == "Requirement"
        return requirement_instrument()
    else
        error("Instrument type not recognized")
    end
end

function requirement_instrument()
    inst = build_instrument(
        FPA_quantum_efficiency = 0.80,
        bench_efficiency = 0.72,
        readout_noise = 120.0,
        dark_current = 10e3u"1/s",
        SSI = 1.0u"nm",
        FWHM = 2.5,
        lower_wavelength = 2040.0,
        upper_wavelength = 2368.0,
        pixel_size_global = 400u"m",
        pixel_size_target = 50u"m"
    )
    return inst
end

function cbe_instrument()
    inst = build_instrument(
        FPA_quantum_efficiency = 0.85,
        bench_efficiency = 0.72,
        readout_noise = 100.0,
        dark_current = 5e3u"1/s",
        SSI = 0.7u"nm",
        FWHM = 1.93,
        lower_wavelength = 2036.0,
        upper_wavelength = 2372.0,
        pixel_size_global = 345u"m",
        pixel_size_target = 34.5u"m"
    )
    return inst
end
