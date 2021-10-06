using CSV
using Plots
using Unitful
using UnitfulRecipes
using Measurements

# We're going to overengineer a bit.

struct Calibration
    wavelength::Vector{typeof(1.0u"nm")}
    counterAvg::Vector{Measurement}
    gradient
    yintercept
end

"""Load calibration data from path string and construct Calibration object."""
function Calibration(path::String)
    headers = ["order", "wavelength", "energy", "min", "max", "avg", "err"] # Custom column headers
    calData = CSV.File(path, header=headers, skipto=2) # skip to first row of data
    wavelength = Quantity.(calData.:wavelength, u"nm")
    counterAvg = calData.:avg .± calData.:err

    # Let's use linear least squares to find the fit
    # Can fit to wavelength = m1 * counts + m2
    # This makes our data kernel:
    kernel = [counterAvg ones(length(counterAvg))]
    # Generalised inverse
    genInv = inv(transpose(kernel) * kernel) * transpose(kernel)
    parameters = genInv * wavelength
    #display(parameters)

    Calibration(wavelength, counterAvg, parameters...)
end

# Expects to be passed an object of type Calibration
@userplot CalibrationPlot
@recipe function f(i::CalibrationPlot)
    cal = i.args[1]
    rawX = cal.counterAvg
    rawY = cal.wavelength

    legend := :bottomright
    xguide := "Counts"
    yguide := "Wavelength"
    grid := true
    minorgrid := true
    seriescolor := :auto

    @series begin
        seriestype := :scatter
        label := "Measurements"
        rawX, rawY
    end

    @series begin
        seriestype := :line
        label := "Fit"
        linestyle := :dash
        y = cal.gradient * rawX .+ cal.yintercept
        rawX, y
    end

end

cal = Calibration("data_calibration_v2.csv")
display(calibrationplot(cal))

λ(count, cal::Calibration=cal) = cal.gradient * count + cal.yintercept
E(λ) = uconvert(u"eV", u"h"* u"c0" / λ)

#------------------------------------------------


data = CSV.File("data_sample_v2.csv")

mutable struct OpticalData
    counter::Vector{Int}
    wavelength::Vector{typeof(1.0u"nm")}
    energy::Vector{typeof((1.0 ± 0.1)u"eV")}
    RawSi::Vector{typeof(1.0u"mV")}
    RawInP::Vector{typeof(1.0u"mV")}
    ScaledSi::Vector{<:Real}
    ScaledInP::Vector{<:Real}
end

"""Load data from path string and construct OpticalData object."""
function OpticalData(path::String)
    headers = ["counter", "wavelength", "energy", "spacer", "withoutInP", "withoutInPScaled", "withInP", "spacer2", "withInPScaled"] # Custom column headers
    data = CSV.File(path, header=headers, skipto=2) # skip to first row of data

    counter = data.:counter
    wavelength = Quantity.(data.:wavelength, u"nm")
    energy = E.(λ.(counter))
    RawSi = Quantity.(data.:withoutInP, u"mV")
    RawInP = Quantity.(data.:withInP, u"mV")
    ScaledSi = RawSi ./ maximum(RawSi)
    ScaledInP = RawInP ./ maximum(RawInP)

    OpticalData(counter, wavelength, energy, RawSi, RawInP, ScaledSi, ScaledInP)
end

data = OpticalData("data_sample_v2.csv")

# Expects to be passed x=energy, y=relative intensity and an optional label or any other plot options
@userplot OpticalDataPlot
@recipe function f(i::OpticalDataPlot)
    x, y = i.args

    legend := :topright
    xguide := "Incident Energy"
    yguide := "Relative Intensity"
    grid := true
    minorgrid := true

    x, y
end

SiPlot = opticaldataplot(data.energy, data.ScaledSi, label="Si")
InPPlot = opticaldataplot(data.energy, data.ScaledInP, label="InP", xlims=(-Inf, 1.5))
display(SiPlot)
display(InPPlot)

@userplot DirectnessCheck
@recipe function f(i::DirectnessCheck)
    energy, intensity = i.args
    invLog = log.(data.RawSi ./ intensity)

    legend := :topleft
    xguide := "Incident Energy"
    yguide := "f(αx)"
    grid := true
    minorgrid := true

    @series begin
        label := "Direct"
        energy, invLog .^ 2
    end

    @series begin
        label := "Indirect"
        energy, sqrt.(invLog)
    end

    @series begin
        # Filter out infinity and NaN rows
        arr = [energy invLog]
        replace!(arr, Inf=>NaN)
        nanrows = any(isnan.(arr), dims=2)
        arr = arr[.!vec(nanrows), :]
        arr[:,2] = arr[:,2] .^ 2

        line = arr[arr[:,1] .> 1.34u"eV", :]
        kernel = [ustrip.(line[:,1]) ones(length(line[:,1]))]
        genInv = inv(transpose(kernel) * kernel) * transpose(kernel)
        parameters = genInv * line[:,2]
        display(parameters)

        fline(x) = parameters[1].val * x + parameters[2].val
        x = LinRange(1.35, 1.385, 2)
        y = fline.(x)

        linestyle := :dash
        label := "Fit"

        x, y
    end
    
end

LinearPlot = directnesscheck(data.energy, data.RawInP, xlims=(1.1, 1.4))
display(LinearPlot)


# -------------------------------------------------