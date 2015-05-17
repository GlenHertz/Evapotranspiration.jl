# Formulas from: 
# http://www.fao.org/docrep/x0490e/x0490e07.htm

module Evapotranspiration

export Location, Angle, Degree, Radian,
       Kelvin, Celcius, Fahrenheit,
       RH, kPa, mm_day, MJ_m2hr, MJ_m2day, Meter,
       ave

abstract Angle
immutable Degree <: Angle
  value::Float64
end
immutable Radian <: Angle
  value::Float64
end
value(a::Angle) = a.value
Base.convert(::Type{Radian}, x::Degree) = π/180*value(x)
Base.promote_rule(::Type{Radian}, ::Type{Degree}) = Radian
Base.full{T<:Angle}(x::T) = x
Base.isnan{T<:Angle}(x::T) = isnan(x.value)
Base.isless{T<:Angle}(x::T, y::T) = isless(x,y)
Base.length{T<:Angle}(x::T) = 1
Base.zero{T<:Angle}(::Type{T}) = T(0.0)
Base.eltype{T<:Angle}(x::T) = eltype(value(x))
Base.getindex{T<:Angle}(x::T, i) = getindex(value(x), i)
Base.abs{T<:Angle}(x::T) = abs(value(x))
Base.sin(x::Radian) = sin(value(x))
Base.sin(x::Degree) = sind(value(x))
Base.cos(x::Radian) = cos(value(x))
Base.cos(x::Degree) = cosd(value(x))
Base.tan(x::Radian) = tan(value(x))
Base.tan(x::Degree) = tand(value(x))
for op in (:+, :-)
  @eval ($op){T<:Angle}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:*,)
  @eval ($op){T<:Radian,N<:Number}(x::N, y::T) = ($op)(x, y.value)
  @eval ($op){T<:Radian,N<:Number}(x::T, y::N) = ($op)(x.value, y)
end

abstract Temperature
immutable Celcius <: Temperature
  value::Float64
end
immutable Kelvin <: Temperature
  value::Float64
end
immutable Fahrenheit <: Temperature
  value::Float64
end
value(a::Temperature) = a.value
Base.convert(::Type{Kelvin},  x::Celcius)    = Kelvin(x.value + 273.15)
#Base.convert(::Type{Celcius}, x::Fahrenheit) = Celcius((x.value-32)/1.8)
Base.convert(::Type{Kelvin},  x::Fahrenheit) = Kelvin(1.8x.value - 457.67)
Kelvin(T::Celcius) = convert(Kelvin, T)
Kelvin(T::Fahrenheit) = convert(Kelvin, T)
Base.promote_rule(::Type{Kelvin}, ::Type{Celcius}) = Kelvin
Base.promote_rule(::Type{Kelvin}, ::Type{Fahrenheit}) = Kelvin
Base.promote_rule(::Type{Kelvin}, ::Type{Float64}) = Float64
Base.full{T<:Temperature}(x::T) = x
Base.isnan{T<:Temperature}(x::T) = isnan(x.value)
Base.isless{T<:Temperature}(x::T, y::T) = isless(x,y)
Base.length{T<:Temperature}(x::T) = 1
Base.zero{T<:Temperature}(::Type{T}) = T(0.0)
Base.eltype{T<:Temperature}(x::T) = eltype(value(x))
Base.getindex{T<:Temperature}(x::T, i) = getindex(value(x), i)
for op in (:+, :-)
  @eval ($op){T<:Temperature}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:*,)
  @eval ($op){T<:Kelvin}(x::T, y::T) = ($op)(x.value, y.value)
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:Temperature}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:Temperature}(x::T, y::Number) = ($op)(x.value, y)
end

immutable RH
  value::Float64
end
value(a::RH) = a.value
for op in (:+, :-)
  @eval ($op){T<:RH}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:RH}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:RH}(x::T, y::Number) = ($op)(x.value, y)
end

immutable MJ_m2day  # MJ/(m^2*day)  (Radiation)
  value::Float64
end
value(a::MJ_m2day) = a.value
for op in (:+, :-)
  @eval ($op){T<:MJ_m2day}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:/,)
  @eval ($op){T<:MJ_m2day}(x::T, y::T) = ($op)(x.value, y.value)
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:MJ_m2day}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:MJ_m2day}(x::T, y::Number) = ($op)(x.value, y)
end
Base.full{T<:MJ_m2day}(x::T) = x
Base.isnan{T<:MJ_m2day}(x::T) = isnan(x.value)
Base.isless{T<:MJ_m2day}(x::T, y::T) = isless(x,y)
Base.length{T<:MJ_m2day}(x::T) = 1
Base.zero{T<:MJ_m2day}(::Type{T}) = T(0.0)
Base.eltype{T<:MJ_m2day}(x::T) = eltype(value(x))
Base.getindex{T<:MJ_m2day}(x::T, i) = getindex(value(x), i)

immutable MJ_m2hr  # MJ/(m^2*hr)  (Radiation)
  value::Float64
end
value(a::MJ_m2hr) = a.value
for op in (:+, :-)
  @eval ($op){T<:MJ_m2hr}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:/,)
  @eval ($op){T<:MJ_m2hr}(x::T, y::T) = ($op)(x.value, y.value)
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:MJ_m2hr}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:MJ_m2hr}(x::T, y::Number) = ($op)(x.value, y)
end
Base.full{T<:MJ_m2hr}(x::T) = x
Base.isnan{T<:MJ_m2hr}(x::T) = isnan(x.value)
Base.isless{T<:MJ_m2hr}(x::T, y::T) = isless(x,y)
Base.length{T<:MJ_m2hr}(x::T) = 1
Base.zero{T<:MJ_m2hr}(::Type{T}) = T(0.0)
Base.eltype{T<:MJ_m2hr}(x::T) = eltype(value(x))
Base.getindex{T<:MJ_m2hr}(x::T, i) = getindex(value(x), i)

immutable mm_day  # mm/day (Evapotranspiration)
  value::Float64
end
value(a::mm_day) = a.value
for op in (:+, :-)
  @eval ($op){T<:mm_day}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:mm_day}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:mm_day}(x::T, y::Number) = ($op)(x.value, y)
end
Base.full{T<:mm_day}(x::T) = x
Base.isnan{T<:mm_day}(x::T) = isnan(x.value)
Base.isless{T<:mm_day}(x::T, y::T) = isless(x,y)
Base.length{T<:mm_day}(x::T) = 1
Base.zero{T<:mm_day}(::Type{T}) = T(0.0)
Base.eltype{T<:mm_day}(x::T) = eltype(value(x))
Base.getindex{T<:mm_day}(x::T, i) = getindex(value(x), i)

immutable kPa  # kPa (1000 Pascals)
  value::Float64
end
value(a::kPa) = a.value
for op in (:+, :-)
  @eval ($op){T<:kPa}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:kPa}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:kPa}(x::T, y::Number) = ($op)(x.value, y)
end
Base.full{T<:kPa}(x::T) = x
Base.isnan{T<:kPa}(x::T) = isnan(x.value)
Base.isless{T<:kPa}(x::T, y::T) = isless(x,y)
Base.length{T<:kPa}(x::T) = 1
Base.zero{T<:kPa}(::Type{T}) = T(0.0)
Base.eltype{T<:kPa}(x::T) = eltype(value(x))
Base.getindex{T<:kPa}(x::T, i) = getindex(value(x), i)

immutable Meter
  value::Float64
end
value(a::Meter) = a.value
for op in (:+, :-)
  @eval ($op){T<:Meter}(x::T, y::T) = T(($op)(x.value, y.value))
end
for op in (:+, :-, :*, :/)
  @eval ($op){T<:Meter}(x::Number, y::T) = ($op)(x, y.value)
  @eval ($op){T<:Meter}(x::T, y::Number) = ($op)(x.value, y)
end
Base.full{T<:Meter}(x::T) = x
Base.isnan{T<:Meter}(x::T) = isnan(x.value)
Base.isless{T<:Meter}(x::T, y::T) = isless(x,y)
Base.length{T<:Meter}(x::T) = 1
Base.zero{T<:Meter}(::Type{T}) = T(0.0)
Base.eltype{T<:Meter}(x::T) = eltype(value(x))
Base.getindex{T<:Meter}(x::T, i) = getindex(value(x), i)


for op in (:*,)
  @eval ($op)(x::kPa, y::RH) = kPa(($op)(x.value, y.value))
  @eval ($op)(x::RH, y::kPa) = kPa(($op)(x.value, y.value))
end



immutable Location
  latitude::Radian  # above (+ve) or below (-ve) the equator
  Lz::Degree # Lz = longitude of the centre of the local time zone [degrees west of Greenwich]
  Lm::Degree # Lm = longitude of the measurement site [degrees west of Greenwich]
  z::Meter   # elevation above sea level
  # Angstrom's values:
  as::Float64 # regression constant 
  bs::Float64
  albedo::Float64
  function Location(lat, Lz, Lm, z=Meter(100), as=0.25, bs=0.5, albedo=0.23)
     new(lat, Lz, Lm, z, as, bs, albedo)
  end
end
function Location(latitude::Degree, Lz::Degree, Lm::Degree, z::Meter=Meter(100), as=0.25, bs=0.5, albedo=0.23) 
  Location(radians(latitude), Lz, Lm, z, as, bs, albedo)
end
latitude(loc::Location) = loc.latitude
Lz(loc::Location) = loc.Lz
Lm(loc::Location) = loc.Lm
Lzm(loc::Location) = value(Lz(loc) - Lm(loc))
elevation(loc::Location) = loc.z
as(loc::Location) = loc.as
bs(loc::Location) = loc.bs
albedo(loc::Location) = loc.albedo
  

using Dates

if false
if !haskey(ENV, "FORECAST_API_KEY")
  error("tests expect a FORECAST_API_KEY to be set in your environment variables")
end

options = {"APIKey" => ENV["FORECAST_API_KEY"],
           "units" => "ca"}
end

ave(x::Number, y::Number) = (x+y)/2
ave{T}(x::T, y::T) = T((x+y)/2)
# Midpoint of hours
function ave(::Type{Hour}, hr1, hr2)
    if hr2 < hr1  # time wrapped around
        hr2 += 24
    end
    (hr1 + hr2)/2
end

# Atmospheric pressure (P)
# P =  atmospheric pressure [kPa],
# z = elevation above sea level [m],
# T = temperature [°C]
P(z::Meter, T::Kelvin) = kPa(101.3((T-0.0065z)/T)^5.26)   # (7)
P(z::Meter, T::Temperature=Celcius(20)) = P(z, Kelvin(T))



## psychrometric constant, gamma, kPa/C
# γ  = psychrometric constant [kPa/°C],
# p  = atmospheric pressure [kPa],
# λ  = latent heat of vaporization, 2.45 [MJ/kg],
# cp = specific heat at constant pressure, 1.013e-3 [MJ/kg/°C],
# e  = ratio molecular weight of water vapour/dry air = 0.622.
γ(z::Meter, T::Kelvin, e=0.622, cp=1.013e-3, λ=2.54) = cp*P(z,T)/(e*λ)  # (8)
γ(p::kPa, e=0.622, cp=1.013e-3, λ=2.54) = cp*p/(e*λ)  # (8)
γ(z::Meter, T::Temperature=Celcius(20), args...) = γ(z, convert(Kelvin, T), args...)

## Relativive humidity, Rh
# ea = actual vapour pressure [kPa],
# e° = saturation vapour pressure at temperature T [kPa],
RH(ea::kPa,T=20) = 100ea/e°(T)   # (10)


## Saturation vapour pressure, e°, in kPa:
# e° = saturation vapour pressure at temperature T [kPa],
# T  = air temperature [°C],
e°(T::Celcius) = kPa(0.6108*exp(17.27T/(T+237.3)))  # (11)

## Mean saturation vapour pressure, Es, in kPa:
# Julia 0.3 bug: would like to use eₛ (instead of Es) but it doesn't work.
Es(Tmin::Celcius, Tmax::Celcius) = ave(e°(Tmax), e°(Tmin))   # (12)

# Slope of saturation vapour pressure curve (Δ) [kPa/C]
Δ(T::Celcius) = 4098*(0.6108*exp(17.27T/(T+237.3)))/(T+237.3)^2 # (13)
Δ(Tmin::Celcius, Tmax::Celcius) = Δ(ave(Tmin, Tmax))

# Actual vapour pressure (Ea) derived from dewpoint temperature
Ea(Tdew::Celcius) = e°(Tdew) # (14)
Ea(Tdew::Temperature) = eᵒ(convert(Celcius, Tdew))

# daily mean RH comes from forecast.io
# This is preferred over Ea(Tdew)
Ea(Tmin::Celcius, Tmax::Celcius, RHmean::RH) = RHmean*ave(e°(Tmax),e°(Tmin)) # (19)
Ea(Tmin::Temperature, Tmax::Temperature, RHmean::RH) = Ea(Celcius(Tmin), Celcius(Tmax), RHmean)

# This is preferred over other methods
Ea(Tmin::Celcius, Tmax::Celcius, RHmin::RH, RHmax::RH) = ave(RHmin*e°(Tmax), RHmax*e°(Tmin)) # (17)

## Vapor pressure deficit, Es - Ea, kPa:
Esa(Tmin::Celcius, Tmax::Celcius; Tdew::Celcius=error("No Tdew")) = Es(Tmin, Tmax) - Ea(Tdew)
Esa(Tmin::Celcius, Tmax::Celcius, RHmean::RH) = Es(Tmin, Tmax) - Ea(Tmin, Tmax, RHmean)
Esa(Tmin::Celcius, Tmax::Celcius, RHmin::RH, RHmax::RH) = Es(Tmin, Tmax) - Ea(Tmin, Tmax, RHmin, RHmax)

######################
# Radiation:
######################

# Extraterrestial radiation

MJ2mm(Ra::MJ_m2day) = MJ_m2day(value(Ra)/2.45)  # (20)  # converts MJ to mm (over whatever time period (eg Day or Hour))
MJ2mm(Ra::MJ_m2hr) = MJ_m2day(value(Ra)/2.45)

# Maximum extraterrestrial radiation for daily periods, Ra
# Ra  = extraterrestrial radiation [MJ/m^2/day],
# ϕ   = latitude [rad] (Equation 22),
# dᵣ  = inverse relative distance Earth-Sun (Equation 23),
# δ   = solar decimation (Equation 24) [rad].
# ωs  = sunset hour angle (Equation 25) [rad],
const Gsc = 0.0820 # Solar constant, Gsc, in MJ/(m^2*min):
function _Ra(ϕ::Radian, dᵣ, δ::Radian, ωs::Radian)
  MJ_m2day(24*60/π * Gsc * dᵣ * (ωs*sin(ϕ)*sin(δ) + cos(ϕ)*cos(δ)*sin(ωs)))  # (21)
end
# latitude = in [deg]
function Ra(::Type{Day}, J, loc::Location)  # (21)
  ϕ = latitude(loc)
  d = δ(J)
  _Ra(ϕ, dᵣ(J), d, ωs(ϕ, d))
end

radians(deg::Degree) = Radian(π/180*value(deg))   # (22)
radians(loc::Location) = radians(latitude(loc))
# dᵣ inverse relative distance Earth-Sun, dᵣ, unitless
dᵣ(J::Day) = 1 + 0.033cos(2π/365*Dates.value(J))  # (23)
# solar decimation, δ
δ(J::Day) = Radian(0.409*sin(2π/365*Dates.value(J)-1.39))  # (24)
# sunset hour angle, ws, in radians (zero is solar noon)
ωs(ϕ::Radian,δ::Radian) = Radian(acos(-tan(ϕ)*tan(δ)))   # (25)
# sun angle at midpoint of period shorter than a day, ws, in radians (zero is solar noon)
# t = standard clock time at the midpoint of the period [hour]. 
#   - For example for a period between 14.00 and 15.00 hours, t = 14.5
# Lz = longitude of the centre of the local time zone [degrees west of Greenwich]
# Lm = longitude of the measurement site [degrees west of Greenwich]
ωs(::Type{Hour}, t, J::Day, loc::Location) = π/12*((t+0.06667*Lzm(loc) + Sc(J))-12)  # (31)
# Sc seasonal correction for solar time [hour]:
function Sc(J::Day)
    b = 2π*(Dates.value(J)-81)/364   # (33)
    0.1645*sin(2b)-0.1255*cos(b)-0.025*sin(b)  # (32)
end



# Max daylight hours
daylighthours(J::Day, loc::Location) = 24/π*ωs(latitude(loc), δ(J))  # (34)
# n = average hours of sun per day
cloudcover(n, J::Day, loc::Location) = 1-n/daylighthours(J,loc)

function duration(::Type{Hour}, hr1, hr2)
    if hr2 < hr1  # time wrapped around
        hr2 += 24
    end
    hr2 - hr1
end
# hr is clock hrs (at midpoint of time period)
function issunset(T::Type{Hour}, hr, J::Day, loc::Location)
    wsun = ωs(T, hr, J, loc)  # 0 rad is noon
    ϕ = latitude(loc)
    d = δ(J)
    # sunset hour angle
    wsunset = ωs(ϕ, d)
    abs(wsun) > abs(wsunset)
end
# Extraterrestrial radiation for hourly periods, Ra, in MJ/(m^2*hour)
function Ra(T::Type{Hour}, J::Day, hr1, hr2, loc::Location) # (28)
    # angle above (below) equator
    ϕ = latitude(loc)
    # solar decimation, d, radians
    d = δ(J)
    # midpiont of hours, eg if 14h and 15h then t=14.5
    hr = ave(T, hr1, hr2)
    # sunset hour angle, in radians
    if issunset(T, hr, J, loc)
        # sun is below horizon:
        return ra = MJ_m2hr(0)
    end
    # The solar time angles at the beginning and end of the period are given by:
    t1 = duration(Hour, hr1, hr2)
    wsun = ωs(T, hr, J, loc)  # sun angle, 0 rad is noon
    w1 = wsun - π*t1/24   # (29)
    w2 = wsun + π*t1/24   # (30)
    ra = 12*60/π*Gsc*dᵣ(J)*((w2-w1)*sin(ϕ)*sin(d)+cos(ϕ)*cos(d)*(sin(w2)-sin(w1)))  # (28)
    MJ_m2hr(ra)
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
function Rs{T}(cloud_cover, J::Day, Ra::T, loc::Location)  # (35)
    N = daylighthours(J, loc)
    sun_pct = (1 - cloud_cover)
    coef = as(loc) + bs(loc)*sun_pct
    rs = T(coef*Ra)
end

# Clear-sky (n ≡ N) solar radiation, Rso
# When as and bs are available:
Rso(Ra,as,bs) = (as+bs)*Ra  # (36)
# When as and bs are not available:
# z = station elevation above sea level [m]
Rso(Ra,z) = (0.75 + 2e-5*z)*Ra  # (37)


# Rs/Rso: Ratio of solar radiation that reaches the earth (compared to a clear sky day)
function Rsso(cloud_cover, J::Day, Ra, loc)
    r1 = Rs(cloud_cover, J, Ra, loc)
    r2 = Rs(0,           J, Ra, loc)
    r1/r2
end

# Net solar radiation
# a = albedo or canopy reflection coefficient
# For example a=0.23 for the hypothetical grass reference crop [dimensionless],
Rns(Rs, loc::Location) = MJ_m2day((1-albedo(loc))Rs)  # (38)

# Net long wave radiation, Rnl [MJ/(m^*day)]
# Rnl = net outgoing longwave radiation [MJ m-2 day-1],
# σ = Stefan-Boltzmann constant [4.903 10-9 MJ K-4 m-2 day-1],
# Tmax = maximum temperature during the 24-hour period [°C],
# Tmin = minimum temperature during the 24-hour period [°C],
# Ea   = actual vapour pressure [kPa],
# Rs/Rso = relative shortwave radiation, 0 <= Rsso <= 1
# Rs measured or calculated. (Equation 35) solar radiation [MJ m-2 day-1],
# Rso calculated (Equation 36 or 37) clear-sky radiation [MJ m-2 day-1].
function Rnl(Rsso, Tmin::Kelvin, Tmax::Kelvin, Ea::kPa)  # (39)
    σ = 4.903e-9
    blackbody = σ * ave(Tmax^4, Tmin^4)
    air_humidity_correction = 0.34 - 0.14*sqrt(value(Ea))
    effect_of_cloudiness = 1.35 * Rsso - 0.351
    MJ_m2day(blackbody * air_humidity_correction * effect_of_cloudiness)
end
function Rnl(Rsso, Tmin::Temperature, Tmax::Temperature, Ea::kPa)
    Rnl(Rsso, convert(Kelvin, Tmin), convert(Kelvin, Tmax), Ea::kPa)
end

# Daily radiation, Rn [MJ/(m^*day)]: 
Rn(Rns, Rnl) = Rns - Rnl   # (40)

function net_radiation(T::Type{Day}, cloud_cover, J::Day, Tmin::Temperature, Tmax::Temperature, RHmin::RH, RHmax::RH, loc::Location)
  ra   = Ra(T, J, loc)
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmin, RHmax)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
end
function net_radiation(T::Type{Day}, cloud_cover, J::Day, Tmin::Temperature, Tmax::Temperature, RHmean::RH, loc::Location)
  ra   = Ra(T, J, loc)
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmean)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
end

# Net radiation [MJ/(m^2*hour)]
function net_radiation(T::Type{Hour}, cloud_cover, J::Day, Tmin, Tmax, RHmin::RH, RHmax::RH, hr1, hr2, loc::Location)
  ra   = Ra(T, J, hr1, hr2, loc) # (28) [MJ/(m^2*hour)]
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmin, RHmax)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
end
function net_radiation(T::Type{Hour}, cloud_cover, J::Day, Tmin, Tmax, RHmean::RH, hr1, hr2, loc::Location)
  ra   = Ra(T, J, hr1, hr2, loc) # (28) [MJ/(m^2*hour)]
  rs   = Rs(cloud_cover, J, ra, loc)
  rsso = Rsso(cloud_cover, J, ra, loc)
  ea   = Ea(Tmin, Tmax, RHmean)
  rnl  = Rnl(rsso, Tmin, Tmax, ea)
  rns  = Rns(rs, loc)
  rn   = Rn(rns, rnl)
end


# Soil heat flux (inaccurate method for over many days)
# G  = soil heat flux [MJ/(m^2*day)],
# cs = soil heat capacity [MJ/(m^3*°C)],
# T1 = air temperature at time 1 [°C],
# T2 = air temperature at time 2 [°C],
# Δt = length of time interval between T1 and T2 [day],
# Δz = effective soil depth [m].
# For a day Δz is between 0.1 and 0.2, 1 month is about 1, 4 months can be 2.0.
x = [1, 30, 4*30]
y = [0.1, 1, 2]
const Δz_slope = (x\y)[1]  # line fit
const Δz_b = y[1]
Δz(::Type{Day}, days) = Δz_slope * days + Δz_b
function soil_heat_flux{Temp<:Union(Celcius,Kelvin)}(T::Type{Day}, Δt, T1::Temp, T2::Temp, Δz=Δz(T,Δt), soil_heat_capacity=2.1)
    0.408*soil_heat_capacity * (T2-T1) * Δz / Δt   # (41)
end

# For monthly periods
# T1 = mean air temperature of previous month
# T2 = mean air temperature of current month
# T3 = mean air temperature of next month
# Δt = length of time interval between T1 and T2 [month],
function soil_heat_flux(T::Type{Month}, Δt, T1::Celcius, T2::Celcius)  # (43) and (44)
    warming = sign(value(T2-T1))
    G = warming*0.14/Δt*(T2-T1)
end

# For hourly or shorter periods
function soil_heat_flux(T::Type{Hour}, hr1, hr2, Rn, latitude, J, Lz, Lm, T1, T2=T1+1)  # (45)
    warming = sign(T2-T1)
    hr = ave(T, hr1, hr2)
    coef = ifelse(issunset(T, hr, latitude, J, Lz, Lm), 0.5, 0.1)
    G = warming*coef*Rn
end

# Windspeed height compensation (for weather stations)
# u2 = wind speed at 2 m above ground surface [m/s],
# uz = measured wind speed at z m above ground surface [m/s],
# z  = height of measurement above ground surface [m].
function u2(windspeed, z=10)
    # speeds are in m/s or kmph
    return windspeed*4.87/log(67.8z-5.42)
end

# ETo = reference evapotranspiration [mm/day],
# Δ   = slope vapour pressure curve [kPa/°C],
# Rn  = net radiation at the crop surface [MJ/(m^2*day)],
# G   = soil heat flux density [MJ/(m^2*day)],
# γ   = psychrometric constant [kPa/°C].
# T   = mean daily air temperature at 2 m height [K],
# u₂  = wind speed at 2 m height [m/s],
# Es  = saturation vapour pressure [kPa],
# Ea  = actual vapour pressure [kPa],
# Esa = Es - Ea = saturation vapour pressure deficit [kPa],
function ETo(Δ, Rn, G, γ, T::Kelvin, u₂, Esa)
   (0.408*Δ*(Rn-G) + 900γ/value(T)*u₂*Esa)/(Δ+γ*(1+0.34u₂))  # (6)
end
function ETo(Td::Type{Day}, Δt, Tmin::Celcius, Tmax::Celcius, RHmin::RH, RHmax::RH, P_kPa::kPa, u10_kph, cloud_cover, J::Day, loc::Location)
    D = Δ(Tmin,Tmax)
    rn = net_radiation(Td, cloud_cover, J, Tmin, Tmax, RHmin, RHmax, loc)
    G = soil_heat_flux(Td, Δt, Tmin, Tmax)
    g = γ(P_kPa)
    T = Kelvin(ave(Tmin,Tmax))
    u₂ = u2(u10_kph/3.6, 10)
    esa = Esa(Tmin, Tmax, RHmin, RHmax)
    ETo(D, rn, G, g, T, u₂, esa)
end
function ETo(Td::Type{Day}, Δt, Tmin::Celcius, Tmax::Celcius, RHmean::RH, P_kPa::kPa, u10_kph, cloud_cover, J::Day, loc::Location)
    D = Δ(Tmin,Tmax)
    rn = net_radiation(Td, cloud_cover, J, Tmin, Tmax, RHmean, loc)
    G = soil_heat_flux(Td, Δt, Tmin, Tmax)
    g = γ(P_kPa)
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
