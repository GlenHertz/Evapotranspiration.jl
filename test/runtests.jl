println("Starting tests")

using Base.Test
using Evapotranspiration
using Base.Dates
using Unitful
using Unitful.DefaultSymbols

@testset "All" begin
@test uconvert(rad, 360°) ≈ 2π*rad

Evap = Evapotranspiration
latitude  = 45.281258°
longitude = -75.575136°

debug=true
Tmax = 19.46°C
Tmin =  5.90°C
Tdew =  6.25°C
RHmean = RH(0.71)
latitude = 45.1°

@test uconvert(K, Tmax) ≈ 292.6K rtol=0.0001

# Number of days through the year, J, days (1 = January 1)
J = 134u"d"

# Lz longitude of the centre of the local time zone
Lz = 75°

# Lm longitude of the measurement site [degrees west of Greenwich]
Lm = 75.76498°

latent_heat_of_vaporaization = 2.54u"MJ/kg" #  MJ/kg



timestamp = 1368594000
d1 = unix2datetime(timestamp)
@test dayofyear(d1) == 135
d2 = Date(2013,2,27)
@test dayofyear(d2) == 58
d3 = Date(2012,12,31)
@test dayofyear(d3) ≈ 366 atol=1 # leap year

# Example 2:
@test Evap.P(1800m, 20°C) ≈ 81.8 atol=0.1
@test Evap.γ(1800m, 20°C) ≈ 0.054 atol=0.01
@test Evap.γ(100m) ≈ 0.0642 atol=0.001

@test Evap.e°(20) ≈ 2.3383 atol=0.0001
# Example 3:
@test Evap.e°(24.5°C) ≈ 3.075 atol=0.001
@test Evap.e°(15.0°C) ≈ 1.705 atol=0.001
@test Evap.Es(15°C, 24.5°C) ≈ 2.39 atol=0.01  # kPa
@test Evap.Es(17.5°C, 29°C) ≈ 3.003 atol=0.001  # kPa

## Saturation vapour pressure slope (Δ) in kPa/C:
@test Evap.Δ(28°C) ≈ 0.2201u"kPa/°C" rtol=0.001 # kPa/C
@test Evap.Δ(26°C, 30°C) ≈ 0.22008u"kPa/°C" atol=0.00001 # kPa/C

## Actual vapour pressure (Ea) derived from dewpoint temperature
@test Evap.Ea(20°C) ≈ 2.338 atol=0.001 # kPa

# daily mean RH comes from forecast.io
# This is preferred over Ea(Tdew)
@test Evap.Ea(20°C, 20°C, RH(0.5)) ≈ 1.169 atol=0.001 # kPa
@test Evap.Ea(18°C, 25°C, ave(RH(0.54),RH(0.82))) ≈ 1.7788 atol=0.0001

# This is preferred over other methods
@test Evap.Ea(18°C, 25°C, RH(0.54), RH(0.82)) ≈ 1.7015 atol=0.0001 # kPa

# Example 5
Tmin = 18°C; RHmax = RH(0.82)
Tmax = 25°C; RHmin = RH(0.54)
@test Evap.e°(Tmin) ≈ 2.064 atol=0.001
@test Evap.e°(Tmax) ≈ 3.168 atol=0.001
@test Evap.Ea(Tmin, Tmax, RHmin, RHmax) ≈ 1.70 atol=0.01
@test Evap.Ea(Tmin, Tmax, ave(RHmin,RHmax)) ≈ 1.78 atol=0.01
# Example 6
esa = Evap.Esa(Tmin, Tmax, RHmin, RHmax)
@test esa ≈ 0.91 atol=0.01

## Vapor pressure deficit, Es - Ea, kPa:
@test Evap.Esa(18°C, 25°C, Tdew=17.8°C) ≈ 0.5777 atol=0.0001  # kPa
@test Evap.Esa(18°C, 25°C, RH(0.68)) ≈ 0.83708 atol=0.00001  # kPa

Radian(2) == Radian(1)

# Example 8:
## Extraterrestrial radiation for daily periods, Ra, in MJ/(m^2*day)
sept3 = Date(2013,9,3)
latitude = Degree(-20) # degrees, 
loc = Location(latitude, Lz, Lm)
ϕ = Evap.latitude(loc)
@test ϕ ≈ -0.35 atol=0.01
J = Day(dayofyear(sept3))
@test J == Day(246)
dᵣ = Evap.dᵣ(J)
@test dᵣ ≈ 0.985 atol=0.001
δ  = Evap.δ(J)
@test δ  ≈ 0.120 atol=0.001
ωs = Evap.ωs(ϕ, δ)
@test ωs ≈ 1.527 atol=0.001
Ra = Evap.Ra(Day, J, loc) # MJ/(m^2*day)
@test Ra ≈ 32.2 atol=0.1 # MJ/(m^2*day)
Ra_mm_day = Evap.MJ2mm(Ra)
@test Ra_mm_day ≈ 13.1 atol=0.1 # mm/day


@test Evap.ave(Hour, 0,1) ≈ 0.5 atol=0.01
@test Evap.ave(Hour, 24, 1) ≈ 24.5 atol=0.01
@test Evap.duration(Hour, 0,1) ≈ 1 atol=0.01
@test Evap.duration(Hour, 24, 1) ≈ 1 atol=0.01
@test Evap.duration(Hour, 24, 4) ≈ 4 atol=0.01
@test Evap.duration(Hour, 24, 0) ≈ 0 atol=0.01
for hr in 0:6
  @test Evap.issunset(Hour, hr, J, loc) == true
end
for hr in 7:17
  @test Evap.issunset(Hour, hr, J, loc) == false
end
for hr in 18:24
  @test Evap.issunset(Hour, hr, J, loc) == true
end


# Example 9:
N = Evap.daylighthours(J, loc)
@test N ≈ 11.7 atol=0.1 # max daylight hours

t = 13.5
ωs = Evap.ωs(Hour, t, J, loc)
@test ωs ≈ 0.385 atol=0.001
# Extraterrestrial radiation for hourly periods, Ra, in mm/time_period)
hr1 = 13
hr2 = hr1+1
J = Day(134)  # dayofyear
ra_24hr = 0
step = 1
for hr1 in 0:step:24
  hr2 = (hr1+step)%24
  ra = Evap.MJ2mm(Evap.Ra(Hour, J, hr1, hr2, loc))
  ra_24hr += ra
end
@test ra_24hr ≈ Evap.MJ2mm(Evap.Ra(Day, J, loc)) atol=0.2

# Example 10:
# The amount of radiation that penetrates the atmosphere, Rs, units same as Ra
latitude = Degree(-(22 + 54/60))
loc = Location(latitude, Lz, Lm)
J = Day(135)    # day of year
Ra = Evap.Ra(Day, J, loc)
@test Ra ≈ 25.1 atol=0.1
N = Evap.daylighthours(J, loc)
@test N ≈ 10.9 atol=0.1
# in May (31 days), 220 hours of sunshine were recorded
n = 220/31
@test n ≈ 7.1 atol=0.1
cloud_cover = Evap.cloudcover(n, J, loc)
@test cloud_cover ≈ 0.348 atol=0.001
z  = 100m
as = 0.25
bs = 0.5
loc = Location(latitude, Lz, Lm, z, as, bs)
Rs = Evap.Rs(cloud_cover, J, Ra, loc)
@test Rs ≈ 14.5 atol=0.1
@test Evap.MJ2mm(Rs) ≈ 5.9 atol=0.1

Tmax = 25.1
Tmin = 19.1
Tdew = 18.1
ra = Evap.Ra(Day, J, loc)
cloud_cover = (1 - 7.1/10.9)
@test Evap.MJ2mm(Evap.Rs(cloud_cover, J, ra, loc)) ≈ 5.90 atol=0.01

# Example 11:
# Rio de Janeiro
latitude = Degree(-(22+54/60))
loc = Location(latitude, Lz, Lm, z, as, bs)
hrs_sunshine_per_month=220
n = hrs_sunshine_per_month/31
N = Evap.daylighthours(J, loc)
cloud_cover =  1 - n/N
@test cloud_cover ≈ 0.349 atol=0.001
Tmin = 19.1°C
Tmax = 25.1°C
Ea = kPa(2.1)  # kPa
Rso = Evap.Rso(Ra,as,bs)
@test Rso ≈ 18.8 atol=0.1
# Rs/Rso: Ratio of solar radiation that reaches the earth (compared to a clear sky day)
Rsso = Evap.Rsso(cloud_cover, J, Ra, loc)
@test Rsso ≈ 0.77 atol=0.01
## Net long wave radiation, Rnl, in mm/day
Rnl = Evap.Rnl(Rsso, Tmin, Tmax, Ea)
@test Rnl ≈ 3.5 atol=0.1

# Example 12:
## Net solar radiation, MJ/(m^2*day)
albedo = 0.23
loc = Location(latitude, Lz, Lm, z, as, bs, albedo)
Rns = Evap.Rns(Rs, loc)
@test Rns ≈ 11.1 atol=0.1
Rn = Evap.Rn(Rns, Rnl)
@test Rn ≈ 7.6 atol=0.1
@test Evap.MJ2mm(Rn) ≈ 3.1 atol=0.1
ea = kPa(2.1)
# Net radiation, Rn, in mm/day: 
net_ra = Evap.net_radiation(Day, cloud_cover, J, Tmin, Tmax, RHmin, RHmax, loc)
@test net_ra ≈ 7.2 atol=0.1

hr1 = 12
hr2 = hr1+1
net_ra_hr = Evap.net_radiation(Hour, cloud_cover, J, Tmin, Tmax, RHmin, RHmax, hr1, hr2, loc)
@test net_ra_hr ≈ -2.363 atol=0.001


## For a day depth is between 0.1 and 0.2, 1 month is about 1, 4 months can be 2.0.
shf = Evap.soil_heat_flux(Day, 30, 14.1°C, 18.8°C)
@test shf ≈ 0.08451 atol=0.00001
@test Evap.soil_heat_flux(Day, 1, Tmin, Tmax) ≈  0.6048 atol=0.0001
@test Evap.soil_heat_flux(Day, 1, Tmax, Tmin) ≈ -0.6048 atol=0.0001
# Example 13:
# Soil is warming
T1 = 14.1°C # March
T2 = 16.1°C # April
T3 = 18.8°C # May
@test Evap.soil_heat_flux(Month, 2, T1, T3) ≈ 0.33 atol=0.01
@test Evap.soil_heat_flux(Month, 1, T1, T2) ≈ 0.28 atol=0.01


# Example 14:
## Windspeed height compensation (for weather stations)
# Find equiv windspeed at 2m
uz = 3.2  # wind speed at 10 m in m/s
z  = 10
u2 = Evap.u2(uz, z)
@test u2 ≈ 2.4 atol=0.01


## Brussels, 6 July located at 50°48'N and at 100 m above sea level:
latitude = Degree(50 + 48/60)
z = 100m
loc = Location(latitude, Lz, Lm, z, as, bs, albedo)
Tmin = 12.3°C
Tmax = 21.5°C
RHmin = RH(0.63)
RHmax = RH(0.84)
RHmean = ave(RHmin,RHmax)
P_kPa = kPa(100.1)
windspeed_10m_kmph = 10
cloud_cover = 0.43
J = Day(187)
eto_day = Evap.ETo(Day, 1, Tmin, Tmax, RHmin, RHmax, P_kPa, windspeed_10m_kmph, cloud_cover, J::Day, loc)
@test eto_day ≈ 3.7 atol=0.01  # mm/day

# Example 17:
# Given the monthly average climatic data of April of Bangkok (Thailand)
# located at 13°44'N and at an elevation of 2 m
latitude = Degree(13+44/60)
z  = 2m  # elevation above sea level
# Lz longitude of the centre of the local time zone
Lz = Degree(100+29/60)
# Lm longitude of the measurement site [degrees west of Greenwich]
Lm = Lz
loc = Location(latitude, Lz, Lm, z, as, bs, albedo)
dt = Date(2015, 4, 15)
J = Day(dayofyear(dt))
Tmin = 25.6°C
Tmax = 34.8°C
u2 = 2
T2 = 30.2°C  # mean monthly average temp, April
T1 = 29.2°C  # mean monthly average temp, March
P_kPa = Evap.P(z, ave(T1,T2))
@test P_kPa ≈ 101.3 atol=0.1
n  = 8.5 # hrs/day of sunshine
cloud_cover = Evap.cloudcover(n, J, loc)
# Need to figure out what to do about RH values: TODO
eto_day = Evap.ETo(Day, 1, Tmin, Tmax, RHmin, RHmax, P_kPa, windspeed_10m_kmph, cloud_cover, J::Day, loc)

# Todo:

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
#	dayofyear = dayofyear(v.time)
#        eto_day = ETo_day(Tmin, Tmax, RHmean, P_kPa, windspeed_10m_kmph, cloud_cover, latitude, dayofyear)
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


println("Done")
end #testset all
