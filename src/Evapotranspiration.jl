# Formulas from: 
# http://www.fao.org/docrep/x0490e/x0490e07.htm

module Evapotranspiration


using Unitful
#using Unitful.DefaultSymbols
import Unitful: Temperature, Pressure, Length, Time, NoUnits

const Temp = Temperature
const Elevation = Length
const DoY = ClockHrs = Time
"λ = latent heat of vaporization, 2.45 [MJ/kg]"
const λ = 2.54u"MJ/kg"
"cp = specific heat at constant pressure, 1.013e-3 [MJ/kg/°C]"
const cp = 1.013e-3u"MJ/kg/°C"
"e′ = ratio molecular weight of water vapour/dry air = 0.622 [unitless]"
const e′ = 0.622
"Solar constant, Gsc, in MJ/(m^2*min)"
const Gsc = 0.0820u"MJ/m^2/minute"

@derived_dimension LatentHeat dimension(u"MJ/kg")
@derived_dimension SpecificHeat dimension(u"MJ/kg/°C")
@derived_dimension HeatFlux dimension(u"W/m^2")
@derived_dimension Speed dimension(u"m/s")

export Location, ave, RH

import Base: *, /, +, -


struct Location
  latitude::DimensionlessQuantity  # above (+ve) or below (-ve) the equator
  Lz::DimensionlessQuantity # Lz = longitude of the centre of the local time zone [+ve west of Greenwich]
  Lm::DimensionlessQuantity # Lm = longitude of the measurement site [degrees west of Greenwich]
  z::Elevation   # elevation above sea level
  # Angstrom's values:
  as::Float64 # regression constant 
  bs::Float64
  albedo::Float64
  function Location(lat::DimensionlessQuantity, Lz::DimensionlessQuantity, Lm::DimensionlessQuantity, z::Elevation=100m, as=0.25, bs=0.5, albedo=0.23)
     new(lat, Lz, Lm, z, as, bs, albedo)
  end
end
latitude(loc::Location) = loc.latitude
Lz(loc::Location) = loc.Lz
Lm(loc::Location) = loc.Lm
Lzm(loc::Location) = uconvert(u"°", Lz(loc) - Lm(loc))
elevation(loc::Location) = loc.z
as(loc::Location) = loc.as
bs(loc::Location) = loc.bs
albedo(loc::Location) = loc.albedo
  

if false
    if !haskey(ENV, "FORECAST_API_KEY")
      error("tests expect a FORECAST_API_KEY to be set in your environment variables")
    end
    
    options = Dict("APIKey" => ENV["FORECAST_API_KEY"],
                   "units" => "ca")
end

# Helper function
"Extract numeric part from quantity in terms of unit"
uval(unit, quantity) = float(uconvert(unit, quantity).val)

"Average Value"
ave(x::T, y::T) where {T} = T((x+y)/2)
ave(x::T1, y::T2) where {T1,T2} = ave(promote(x, y)...)
# Midpoint of hours
function ave_hrs(t1::ClockHrs, t2::ClockHrs)
    hr1 = uconvert(u"hr", t1)
    hr2 = uconvert(u"hr", t2)
    if hr2 < hr1  # time wrapped around
        hr2 += 24u"hr"
    end
    (hr1 + hr2)/2
end


# Atmospheric pressure (P)
# P = atmospheric pressure [kPa],
# z = elevation above sea level [m],
# T = temperature [K]
"P = Atmospheric pressure, z = elevation, T = temperature, [kPa]"
function P(z::Elevation, T::Temp)
    z = uconvert(u"m", z)
    T = uconvert(u"K", T)
    101.3u"kPa" * ((T - 0.0065u"K/m"*z)/T)^5.26   # (7)
end

## psychrometric constant, gamma, kPa/°C   # (8)
# γ′ = psy = psychrometric constant [kPa/°C],
# p  = atmospheric pressure [kPa],
# λ  = latent heat of vaporization, 2.45 [MJ/kg],
# cp = specific heat at constant pressure, 1.013e-3 [MJ/kg/°C],
# e′  = ratio molecular weight of water vapour/dry air = 0.622.
"γ′ = psychrometric constant, gamma, [kPa/°C]"
function γ′(p::Pressure)
    uconvert(u"kPa/°C", cp*p/(e′*λ))  # (8)
end
"γ′ = psychrometric constant, gamma, [kPa/°C]"
γ′(z::Elevation, t::Temp) = γ′(P(z, t))

## Relativive humidity, Rh
# ea = actual vapour pressure [kPa],
# e° = saturation vapour pressure at temperature t [kPa],
"Relative humidity [unitless]"
function RH(ea::Pressure, t::Temp=20°C)
    ea = uconvert(u"kPa", ea)
    rh = 100ea/e°(t)   # (10)
    uconvert(Unitful.NoUnits, rh)
end


## Saturation vapour pressure, e°, in kPa:
# e° = saturation vapour pressure at temperature t [kPa],
# t  = air temperature [°C],
"e° = Saturation vapour pressure [kPa]"
function e°(t::Temp)
    t = uconvert(u"°C", t)
    0.6108u"kPa" * exp(17.27t/(t+237.3))      # (11)
end

## Mean saturation vapour pressure, Es, in kPa:
# Julia 0.3 bug: would like to use eₛ (instead of Es) but it doesn't work.
"Mean saturation vapour pressure, Es [kPa]"
function Es(Tmin::Temp, Tmax::Temp) 
    p = ave(e°(Tmax), e°(Tmin))   # (12)
    uconvert(u"kPa", p)
end

# Slope of saturation vapour pressure curve (Δ) [kPa/°C]
# T  = mean air temperature [°C]
"Slope of saturation vapour pressure curve (Δ) [kPa/°C]"
function Δ(T::Temp) 
    T = uconvert(u"°C", T)
    c1 = 4098u"°C"
    t1 = T+237.3u"°C"
    c1*(e°(T))/t1^2 # (13)
end
Δ(Tmin::Temp, Tmax::Temp) = Δ(ave(Tmin, Tmax))

"Actual vapour pressure (Ea) derived from dewpoint temperature [kPa]"
Ea(Tdew::Temp) = e°(Tdew) # (14)

# daily mean RH comes from forecast.io
# This is preferred over Ea(Tdew)
"Actual vapour pressure (Ea) derived from dewpoint temperature [kPa]"
function Ea(Tmin::Temp, Tmax::Temp, RHmean::Real)
    RHmean*ave(e°(Tmax),e°(Tmin)) # (19)
end

# This is preferred over other methods
"Actual vapour pressure (Ea) derived from dewpoint temperature [kPa]"
function Ea(Tmin::Temp, Tmax::Temp, RHmin::Real, RHmax::Real)
    ave(e°(Tmin)*RHmax, e°(Tmax)*RHmin) # (17)
end

"Vapor pressure deficit, Es - Ea [kPa]"
function Esa(Tmin::Temp, Tmax::Temp; Tdew::Temp=error("No Tdew"))
    Es(Tmin, Tmax) - Ea(Tdew)
end
"Vapor pressure deficit, Es - Ea [kPa]"
function Esa(Tmin::Temp, Tmax::Temp, RHmean)
    Es(Tmin, Tmax) - Ea(Tmin, Tmax, RHmean)
end
"Vapor pressure deficit, Es - Ea [kPa]"
function Esa(Tmin::Temp, Tmax::Temp, RHmin, RHmax)
    Es(Tmin, Tmax) - Ea(Tmin, Tmax, RHmin, RHmax)
end

######################
# Radiation:
######################

# Extraterrestial radiation
FIXME:
MJ2mm(Ra::MJ_m2day) = MJ_m2day(value(Ra)/2.45)  # (20)  # converts MJ to mm (over whatever time period (eg Day or Hour))
MJ2mm(Ra::MJ_m2hr) = MJ_m2day(value(Ra)/2.45)

# Maximum extraterrestrial radiation for daily periods, Ra
# Ra  = extraterrestrial radiation [MJ/m^2/day],
# ϕ   = latitude [rad] (Equation 22),
# dᵣ  = inverse relative distance Earth-Sun (Equation 23),
# δ   = solar decimation (Equation 24) [rad].
# ωs  = sunset hour angle (Equation 25) [rad],
"Maximum extraterrestrial radiation for daily periods, Ra [MJ/m^2/day]"
function Ra(ϕ::DimensionlessQuantity, dᵣ::Real, δ::DimensionlessQuantity, ωs::DimensionlessQuantity)
  uconvert(u"MJ/m^2/d", (1/π * Gsc * dᵣ * (ωs*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ωs))))  # (21)
end

"Maximum extraterrestrial radiation for daily periods, Ra [MJ/m^2/day]"
function Ra(J::DoY, loc::Location)  # (21)
  ϕ = latitude(loc)
  Ra(ϕ, dᵣ(J), δ(J), ωs(ϕ, d))
end

"dᵣ inverse relative distance Earth-Sun, dᵣ, unitless"
dᵣ(J::DoY) = 1 + 0.033cos(2π/365*uval(u"d", J))  # (23)
function dᵣ(J::DoY)
    J = uconvert(u"d", J)
    r = 2π*u"rad" / 365u"d" * J
    1 + 0.033cos(r)  # (23)
end
"solar decimation, δ, [rad]"
function δ(J::DoY)
    J = uconvert(u"d", J)
    α = 2π * J/365u"d"
    (0.409*sin(α)-1.39)u"rad"  # (24)
end
"sunset hour angle, ws, in radians (zero is solar noon) [rad]"
function ωs(ϕ::DimensionlessQuantity, δ::DimensionlessQuantity)
    (acos(-tan(ϕ)*tan(δ)))u"rad"   # (25)
end
# sun angle at midpoint of period shorter than a day, ws, in radians (zero is solar noon)
# t = standard clock time at the midpoint of the period [hour]. 
#   - For example for a period between 14.00 and 15.00 hours, t = 14.5
# Lz = longitude of the centre of the local time zone [degrees west of Greenwich]
# Lm = longitude of the measurement site [degrees west of Greenwich]
"sun angle at midpoint of period shorter than a day, ws, (zero is solar noon) [rad]"
function ωs(t::ClockHrs, J::DoY, loc::Location) 
    t = uconvert(u"hr", t)
    J = uconvert(u"d", J)
    c1 = 0.06667u"hr"
    π*u"rad"/12u"hr" * ((t + c1*Lzm(loc) + Sc(J)) - 12u"hr")  # (31)
end

"Seasonal correction for solar time, Sc [hour]"
function Sc(J::DoY)
    J = uconvert(u"d", J)
    b = 2π*u"rad" * (J - 81u"d")/364u"d"   # (33)
    (0.1645*sin(2b) - 0.1255*cos(b) - 0.025*sin(b))u"hr"  # (32)
end

"Max daylight hours"
function daylighthours(J::DoY, loc::Location)
    24u"hr"/π*ωs(latitude(loc), δ(J))  # (34)
end

"Cloud cover percent with n = hours of sun in a day, J = day of year [%]"
function cloudcover(n::ClockHrs, J::DoY, loc::Location)
    1 - n/daylighthours(J, loc)
end

"Number of hours between two clock hours"
function hrs_duration(t1::ClockHrs, t2::ClockHrs)
    hr1 = convert(u"hr", t1)
    hr2 = convert(u"hr", t1)
    if hr2 < hr1  # time wrapped around
        hr2 += 24u"hr"
    end
    hr2 - hr1
end

"Check if clock hours is sunset (at midpoint of time period), J = day of year"
function issunset(hr::ClockHrs, J::DoY, loc::Location)
    wsun = ωs(hr, J, loc)  # 0 rad is noon
    ϕ = latitude(loc)
    d = δ(J)
    # sunset hour angle
    wsunset = ωs(ϕ, d)
    abs(wsun) > abs(wsunset)
end

"Percentage of day for duration of two times"
function day_pct(t1::ClockHrs, t2::ClockHrs)
    convert(NoUnits, hrs_duration(t1, t2)/24u"hr")
end

"Extraterrestrial radiation for hourly periods, Ra [MJ/m^2/hour]"
function Ra(J::DoY, t1::ClockHrs, t2::ClockHrs, loc::Location) # (28)
    # angle above (below) equator
    ϕ = latitude(loc)
    # solar decimation, d, radians
    d = δ(J)
    # sunset hour angle, in radians
    tmid = ave_hrs(t1, t2)
    if issunset(tmid, J, loc)
        # sun is below horizon:
        return 0.0u"MJ/m^2/hr"
    end
    # The solar time angles at the beginning and end of the period are given by:
    w1 = ωs(t1, J, loc)  # sun angle, 0 rad is noon
    w2 = ωs(t2, J, loc)  # sun angle, 0 rad is noon
    ra = Gsc*dᵣ(J)*12/π*((w2-w1)*sin(ϕ)*sin(d)+cos(ϕ)*cos(d)*(sin(w2)-sin(w1)))  # (28)
    uconvert(u"MJ/m^2/hr", ra)
end
    
### Solar or Shortwave Radiation, Rs

# The amount of radiation that penetrates the atmosphere, Rs, units same as Ra
# Rs solar or shortwave radiation 
# cloud_cover = percentage of cloud coverage, 0.5 is half covered
# Ra = extraterrestrial radiation (same units as Rs)
# as = regression constant, expressing the fraction of extraterrestrial radiation reaching the earth on overcast days (n = 0),
# as+bs = fraction of extraterrestrial radiation reaching the earth on clear days (n ≡ N (clear sky)).
# n = hours of sun per day
# N = max hours of sun per day
"The amount of radiation that penetrates the atmosphere, Rs, units same as Ra"
function Rs(cloud_cover::Real, J::DoY, Ra::RaT, loc::Location) where {RaT} # (35)
    #?? N = daylighthours(J, loc)
    sun_pct = (1 - cloud_cover)
    coef = as(loc) + bs(loc)*sun_pct
    RaT(coef*Ra)
end

# Clear-sky (n ≡ N) solar radiation, Rso
# When as and bs are available:
"Clear-sky solar radiation (n ≡ N), Rso"
function Rso(Ra::HeatFlux, as::Real, bs::Real)
    (as+bs)*Ra  # (36)
end

# When as and bs are not available:
# z = station elevation above sea level [m]
"Clear-sky solar radiation (n ≡ N) at elevation z, Rso"
function Rso(Ra::HeatFlux, z::Elevation)
    (0.75 + z/50000u"m")*Ra  # (37)
end


"Rs/Rso: Ratio of solar radiation that reaches the earth (compared to a clear sky day)"
function Rsso(cloud_cover::Real, J::DoY, Ra::HeatFlux, loc::Location)
    r1 = Rs(cloud_cover, J, Ra, loc)
    r2 = Rs(0,           J, Ra, loc)
    r1/r2
end

# Net solar radiation
# a = albedo or canopy reflection coefficient
# For example a=0.23 for the hypothetical grass reference crop [dimensionless],
"Net short wave solar radiation, which is a balance between the incoming and reflected [MJ/m^2/d]"
function Rns(Rs::HeatFlux, loc::Location)
    ((1-albedo(loc))*Rs)  # (38)
end

# Net long wave radiation, Rnl [MJ/(m^*day)]
# Rnl = net outgoing longwave radiation [MJ m-2 day-1],
# σ = Stefan-Boltzmann constant [4.903 10-9 MJ K-4 m-2 day-1],
# Tmax = maximum temperature during the 24-hour period [°C],
# Tmin = minimum temperature during the 24-hour period [°C],
# Ea   = actual vapour pressure [kPa],
# Rs/Rso = relative shortwave radiation, 0 <= Rsso <= 1
# Rs measured or calculated. (Equation 35) solar radiation [MJ m-2 day-1],
# Rso calculated (Equation 36 or 37) clear-sky radiation [MJ m-2 day-1].
"Net outgoing long wave radiation, Rnl [MJ/(m^2*day)]"
function Rnl(Rsso::Real, Tmin::Temp, Tmax::Temp, Ea::Pressure)  # (39)
    σ = 4.903e-9u"MJ/K^4/m^2/d"
    blackbody = σ * ave(u"K"(Tmax)^4, u"K"(Tmin)^4)
    air_humidity_correction = 0.34 - 0.14u"1/kPa"*Ea
    effect_of_cloudiness = 1.35 * Rsso - 0.351
    blackbody * air_humidity_correction * effect_of_cloudiness
end

"Daily radiation, Rn [MJ/(m^2*day)]" 
function Rn(Rns::HeatFlux, Rnl::HeatFlux)
    Rns - Rnl   # (40)
end

"Net radiation per day [MJ/m^2/day]"
function net_radiation(cloud_cover::Real, J::DoY, Tmin::Temp, Tmax::Temp, RHmin::Real, RHmax::Real, loc::Location)
  ra   = Ra(J, loc)
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmin, RHmax)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
  uconvert(u"MJ/m^2/d", rn)
end
"Net radiation per day [MJ/m^2/day]"
function net_radiation(cloud_cover::Real, J::DoY, Tmin::Temp, Tmax::Temp, RHmean::Real, loc::Location)
  ra   = Ra(J, loc)
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmean)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
  uconvert(u"MJ/m^2/d", rn)
end

"Net radiation [MJ/(m^2*hour)]"
function net_radiation(cloud_cover::Real, J::DoY, Tmin::Temp, Tmax::Temp, RHmin::Real, RHmax::Real, hr1::ClockHrs, hr2::ClockHrs, loc::Location)
  ra   = Ra(J, hr1, hr2, loc) # (28) [MJ/(m^2*hour)]
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmin, RHmax)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
  uconvert(u"MJ/m^2/hr", rn)
end
"Net radiation [MJ/(m^2*hour)]"
function net_radiation(cloud_cover::Real, J::DoY, Tmin::Temp, Tmax::Temp, RHmean::Real, hr1::ClockHrs, hr2::ClockHrs, loc::Location)
  ra   = Ra(J, hr1, hr2, loc) # (28) [MJ/(m^2*hour)]
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmean)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
  uconvert(u"MJ/m^2/hr", rn)
end

# Soil heat flux (inaccurate method for over many days)
# G  = soil heat flux [MJ/(m^2*day)],
# cs = soil heat capacity [MJ/(m^3*°C)],
# T1 = air temperature at time 1 [°C],
# T2 = air temperature at time 2 [°C],
# Δt = length of time interval between T1 and T2 [day],
# Δz = effective soil depth [m].
# For a day Δz is between 0.1 and 0.2, 1 month is about 1, 4 months can be 2.0.
#x = [1, 30, 4*30] * u"d"
#y = [0.1, 1, 2] * u"m"
"Soil heat flux [MJ/m^2/d]"
function soil_heat_flux(Δt::Time, T1::Temp, T2::Temp)
    if Δt <= 1u"d"
        soil_heat_flux_hr(Δt, T1, T2)
    elseif Δt <= 10u"d"
        0.0u"MJ/m^2/d"
    elseif Δt <= 30u"d" 
        soil_heat_flux_days(Δt, T1, T2)
    else
        soil_heat_flux_months(Δt, T1, T2)
    end
end
"Effective soil depth for heat flux [m]"
function Δz(t::Time)
    days = uconvert(u"d", t)
    if days <= 1u"d"
        return 0.1u"m"
    elseif days <= 30u"d"
        return (1u"m"-0.1u"m")/30u"d" * days + 0.1u"m"
    elseif days <= 4*30u"d"
        return (1u"m")/(3*30u"d") * days + 1u"m"
    else
        return 2.0u"m"
    end
end
"Soil heat flux calculated over days [MJ/m^2/d]"
function soil_heat_flux_days(Δt::Time, T1::Temp, T2::Temp, Δz=Δz(T,Δt))
    soil_heat_capacity=2.1u"MJ/m^3/°C"  # depends on soil composition
    t1 = uconvert(u"°C", T1)
    t2 = uconvert(u"°C", T2)
    Δt = uconvert(u"d",  Δt)
    G = soil_heat_capacity * (t2-t1) * Δz(Δt) / Δt   # (41)
    uconvert(u"MJ/m^2/d", G)
end

# For monthly periods
# T1 = mean air temperature of previous month
# T2 = mean air temperature of current month
# T3 = mean air temperature of next month
# Δt = length of time interval between T1 and T2 [month],
"Soil heat flux calculated over months [MJ/m^2/d]"
function soil_heat_flux_months(Δt::Time, T1::Temp, T2::Temp)  # (43) and (44)
    T1 = uconvert(u"°C", T1)
    T2 = uconvert(u"°C", T2)
    Δt = uconvert(u"d",  Δt)
    warming = sign(T2-T1)
    G = warming*0.14u"MJ/m^2/°C"/Δt*(T2-T1)
    uconvert(u"MJ/m^2/d", G)
end

# For hourly or shorter periods
function soil_heat_flux_hours(hr1::ClockHrs, hr2::ClockHrs, Rn::HeatFlux, J::DoY, loc::Location, T1::Temp, T2::Temp)  # (45)
    warming = sign(T2-T1)
    hr = ave_hrs(hr1, hr2)
    coef = ifelse(issunset(hr, J, loc), 0.5, 0.1)
    G = warming*coef*Rn                           # (45) (46)
    uconvert(u"MJ/m^2/d", G)
end

# Windspeed height compensation (for weather stations)
# u2 = wind speed at 2 m above ground surface [m/s],
# uz = measured wind speed at z m above ground surface [m/s],
# z  = height of measurement above ground surface [m].
function u2(windspeed::Speed, z::Elevation=10u"m")
    # speeds are in m/s or kmph
    c1 = (67.8z-5.42u"m")/1u"m"
    windspeed*4.87/log(c1)                    # (47)
end

# ETo = reference evapotranspiration [mm/day],
# Δ   = slope vapour pressure curve [kPa/°C],
# Rn  = net radiation at the crop surface [MJ/(m^2*day)],
# G   = soil heat flux density [MJ/(m^2*day)],
# γ′   = psychrometric constant [kPa/°C].
# T   = mean daily air temperature at 2 m height [K],
# u₂  = wind speed at 2 m height [m/s],
# Es  = saturation vapour pressure [kPa],
# Ea  = actual vapour pressure [kPa],
# Esa = Es - Ea = saturation vapour pressure deficit [kPa],
function ETo(Δ, Rn, G, γ′, T::typeof(K), u₂, Esa)
   (0.408*Δ*(Rn-G) + 900γ/value(T)*u₂*Esa)/(Δ+γ′*(1+0.34u₂))  # (6)
end
function ETo(Td::Type{Day}, Δt, Tmin::typeof(°C), Tmax::typeof(°C), RHmin, RHmax, P_kPa::typeof(kPa), u10_kph, cloud_cover, J::DoY, loc::Location)
    D = Δ(Tmin,Tmax)
    rn = net_radiation(Td, cloud_cover, J, Tmin, Tmax, RHmin, RHmax, loc)
    G = soil_heat_flux(Td, Δt, Tmin, Tmax)
    g = γ′(P_kPa)
    T = Kelvin(ave(Tmin,Tmax))
    u₂ = u2(u10_kph/3.6, 10)
    esa = Esa(Tmin, Tmax, RHmin, RHmax)
    ETo(D, rn, G, g, T, u₂, esa)
end
function ETo(Td::Type{Day}, Δt, Tmin::typeof(°C), Tmax::typeof(°C), RHmean, P_kPa::typeof(kPa), u10_kph, cloud_cover, J::DoY, loc::Location)
    D = Δ(Tmin,Tmax)
    rn = net_radiation(Td, cloud_cover, J, Tmin, Tmax, RHmean, loc)
    G = soil_heat_flux(Td, Δt, Tmin, Tmax)
    g = γ′(P_kPa)
    T = Kelvin(ave(Tmin,Tmax))
    u₂ = u2(u10_kph/3.6, 10)
    esa = Esa(Tmin, Tmax, RHmean)
    ETo(D, rn, G, g, T, u₂, esa)
end

saturation_vapor_pressure(ave_temp_celcius) = 10^(8.07131-1730.63/(ave_temp_celcius+233.426))

function saturation_vapor_pressure_slope(ave_temp_celcius, dt=0.1)
    p1 = saturation_vapor_pressure(ave_temp_celcius)
    p2 = saturation_vapor_pressure(ave_temp_celcius+dt*5/9)
    return (p2-p1)/dt
end

temperature_coef(ave_temp_celcius, scaling=0.03567) = scaling*5.67e-8*(273+ave_temp_celcius)^4

function evaporation_mm_per_day(wind_speed_kmph, ave_temp_celcius, dew_point_celcius)
    mi_per_km = 1/1.61
    wind_speed_mi_per_day = wind_speed_kmph * 24 * mi_per_km
    ea = saturation_vapor_pressure(ave_temp_celcius)
    ed = saturation_vapor_pressure(dew_point_celcius)
    0.35*(ea-ed)*(1+0.0098*wind_speed_mi_per_day)
end

function degree_days_mm_per_day(ave_temp_celcius, dew_point_celcius, solar_radiation, sunshine_pct)
    mi_per_km = 1/1.61
    Ra = solar_radiation
    Re = 0.15 # reflexivity of forrest
    T = ave_temp_celcius
    Tcoef = 0.03567*5.67e-8*(273+T)^4
    ed = saturation_vapor_pressure(dew_point_celcius)
    S = sunshine_pct
    (Ra*(1-Re)*(0.18+0.55*S))-(Tcoef*(0.56-0.092*ed^0.5)*(0.1+0.9*S))
end

function potential_evaporation(ave_temp_celcius, dew_point_celcius, wind_speed_kmph, solar_radiation, sunshine_pct)
    Slope = saturation_vapor_pressure_slope(ave_temp_celcius)
    E = evaporation_mm_per_day(wind_speed_kmph, ave_temp_celcius, dew_point_celcius)
    H = degree_days_mm_per_day(ave_temp_celcius, dew_point_celcius, solar_radiation, sunshine_pct)
    ((Slope*H)+(0.27*E))/(Slope+0.27)
end

end # module

#function print_evaporation(forecast_data) {
#    vec = forecast_data
#    totalDewPoint = 0
#    totalRH = 0
#    numDewPoints = 0
#    for (i=0; i < vec.length; i++) {
#        v = vec[i]
#	v.date = Date(v.time*1000)
#	v.dom = v.date.getDate()
#	v.hr = v.date.getHours()
#        ave_temp_c = NaN
#        if ("temperature" in v) {
#            ave_temp_c = v.temperature
#        } else if (("temperatureMax" in v) && ("temperatureMin" in v)) {
#            ave_temp_c = v.temperatureMin + 0.75*(v.temperatureMax - v.temperatureMin)
#        }
#        dew_temp_c = v.dewPoint
#        solar_radiation = 6
#        sunshine = 1-v.cloudCover; # account for night having no sun.
#        wind_speed_kmph = v.windSpeed
#        E = evaporation_mm_per_day(wind_speed_kmph, ave_temp_c, dew_temp_c)
#        H = degree_days_mm_per_day(ave_temp_c, dew_temp_c, 6, sunshine)
#                U = potential_evaporation(ave_temp_c, dew_temp_c, wind_speed_kmph, solar_radiation, sunshine)
#        precipIntensity = NaN
#        if ("percipIntensityMax" in v) {
#            precipIntensity = v.percipIntensityMax
#        }
#        if ("precipIntensity" in v) {
#            precipIntensity = v.precipIntensity
#        }
#        percip_mm_day = v.precipIntensity; # * vec[i].precipProbability
#        if (v.dom == 18) {
#		console.log(v.date+': '+v.humidity+', '+v.dewPoint)
#		totalDewPoint += v.dewPoint
#		totalRH += v.humidity
#		numDewPoints += 1
#	}
#        if (typeof percip_mm_day == 'undefined') {
#            percip_mm_day = 0
#        }
#        AM = percip_mm_day - U
#        t = Date(1000*vec[i].time)
#        #console.log(vec[i].time+": "+t + ": E="+ E.toFixed(2) + ", H=" + H.toFixed(2)+", ETo="+U.toFixed(2)+", percip="+percip_mm_day.toFixed(2)+", AM="+AM.toFixed(2))
#
#	# New method:
#	Tmin = v.temperatureMin
#	Tmax = v.temperatureMax
#	Tdew = v.dewPoint
#	RHmean = v.humidity
#	P_kPa = v.pressure
#	# TODO: scale down wind because of shelter from trees, maybe 0.25
#	windspeed_10m_kmph = v.windSpeed
#	cloud_cover = v.cloudCover
#	J = J(v.time)
#        eto_day = ETo_day(Tmin, Tmax, RHmean, P_kPa, windspeed_10m_kmph, cloud_cover, latitude, J)
#	eto_diff = U - eto_day
#        #console.log("ETo day="+eto_day.toFixed(2)+" [mm/day]   diff="+eto_diff)
#    }
#    console.log("Total dew: "+ totalDewPoint+" RH: "+totalRH)
#    aveDewPoint = totalDewPoint/numDewPoints
#    aveRH = totalRH/numDewPoints
#    console.log("Ave: Dew Point: "+ aveDewPoint+", RH: "+ aveRH)
#}
#
#
#function print_daily_forecast(data) {
#  daily = data.daily.data
#  for (i=0; i<daily.length; i++) {
#      di = daily[i]
#      console.log(Date(1000*di.time) + ":")
#      console.log("  - sunrise: " + Date(1000*di.sunriseTime))
#      console.log("  - sunset:  " + Date(1000*di.sunsetTime))
#      console.log("  - TmaxTime: " + Date(1000*di.temperatureMaxTime))
#      console.log("  - TminTime:  " + Date(1000*di.temperatureMinTime))
#      console.log(di)
#  }
#}
#
#forecast.get(latitude, longitude, options, function (err, res, data) {
#  if (err) throw err
#  #console.log('res: ' + util.inspect(res))
#  #console.log('data: ' + util.inspect(data, true, 5))
#  #console.log('daily data: ' + util.inspect(data, true, 5))
#  #print_daily_forecast(data)
#  print_evaporation(data.daily.data)
#  print_evaporation(data.hourly.data)
#})
#
#
