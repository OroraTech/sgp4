package sgp4

import (
	"fmt"
	"math"
	"slices"
	"time"
)

const MinElevationForPass = 10.0 // degrees

type Eci struct {
	DateTime time.Time
	Position Vector
	Velocity Vector
}

type Vector struct {
	X, Y, Z float64
}

// ToGeodetic converts ECI coordinates to geodetic coordinates (lat, lon, alt)
// ToGeodetic converts ECI coordinates to geodetic coordinates (lat, lon, alt)
func (eci *Eci) ToGeodetic() (lat, lon, alt float64) {
	currentRe := reSGP4 // Use SGP4-aligned constants
	currentF := fSGP4
	e2 := currentF * (2.0 - currentF)

	gmst := eci.GreenwichSiderealTime()
	x := eci.Position.X
	y := eci.Position.Y
	z := eci.Position.Z

	lon = math.Atan2(y, x) - gmst
	lon = wrapLongitude(lon)
	r := math.Sqrt(x*x + y*y)
	lat = AcTan(z, r)

	const maxIter = 10
	const tol = 1e-10
	var oldLat float64
	var c_iter float64

	for range maxIter {
		oldLat = lat
		sinLat := math.Sin(lat)
		if math.Abs(1.0-e2*sinLat*sinLat) < 1e-14 { // Avoid division by zero or sqrt of negative
			c_iter = 1.0 / math.Sqrt(1e-14)
		} else {
			c_iter = 1.0 / math.Sqrt(1.0-e2*sinLat*sinLat)
		}
		lat = AcTan(z+currentRe*c_iter*e2*sinLat, r)
		if math.Abs(lat-oldLat) < tol {
			break
		}
	}

	sinLat := math.Sin(lat)
	cosLat := math.Cos(lat)
	var N_val float64
	if math.Abs(1.0-e2*sinLat*sinLat) < 1e-14 {
		N_val = currentRe * (1.0 / math.Sqrt(1e-14))
	} else {
		N_val = currentRe * (1.0 / math.Sqrt(1.0-e2*sinLat*sinLat))
	}

	if math.Abs(cosLat) < 1e-10 {
		alt = math.Abs(z) - currentRe*math.Sqrt(1.0-e2)
	} else {
		alt = r/cosLat - N_val
	}

	lat = lat * rad2deg
	lon = lon * rad2deg
	return lat, lon, alt
}

// FindPosition propagates the TLE to the given time offset (tsince) in minutes
// FindPosition propagates the TLE to the given time offset (tsince) in minutes
func (tle *TLE) FindPosition(tsince float64) (Eci, error) {
	elems, err := tle.Initialize()
	if err != nil {
		return Eci{}, fmt.Errorf("SGP4 propagation error during initialization: %w", err)
	}

	// Retrieve initialized SGP4 elements and constants
	// These are mean elements at epoch, some are already perturbed by J2 (a, n)
	// Others are TLE raw values (ecc, m, omega, raan, incl)
	// And derived SGP4 constants (c1, c4, xmdot, etc.)

	// Secular effects (similar to SGP4::FindPositionSGP4 before CalculateFinalPositionVelocity)
	xmdf := elems.m + elems.xmdot*tsince
	omgadf := elems.omega + elems.omgdot*tsince
	xnoddf := elems.raan + elems.xnodot*tsince

	omega := omgadf // current argument of perigee
	xmp := xmdf     // current mean anomaly

	tsq := tsince * tsince
	xnode := xnoddf + elems.xnodcf*tsq // current RAAN
	tempa := 1.0 - elems.c1*tsince
	tempe := elems.bstar * elems.c4 * tsince
	templ := elems.t2cof * tsq

	if !elems.isSimpleModel {
		delomg := elems.omgcof * tsince
		delm_term := 0.0
		if elems.eta != 0.0 { // Avoid division by zero if eta is zero (circular orbit or specific init issue)
			delm_term = elems.xmcof * (math.Pow(1.0+elems.eta*math.Cos(xmdf), 3.0) - elems.delmo)
		}

		temp := delomg + delm_term
		xmp += temp
		omega -= temp

		tcube := tsq * tsince
		tfour := tsince * tcube
		tempa = tempa - elems.d2*tsq - elems.d3*tcube - elems.d4*tfour
		tempe += elems.bstar * elems.c5 * (math.Sin(xmp) - elems.sinmo)
		templ += elems.t3cof*tcube + tfour*(elems.t4cof+tsince*elems.t5cof)
	}

	a := elems.a * tempa * tempa              // current semi-major axis (ER)
	e := elems.ecc - tempe                    // current eccentricity
	xl := xmp + omega + xnode + elems.n*templ // current mean longitude (M + omega + Omega)

	// Ensure eccentricity is within sane bounds
	if e <= -0.001 { // Check current eccentricity 'e'
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonEccentricityTooLow, Value: e}
	} else if e < 1.0e-6 {
		e = 1.0e-6
	} else if e > (1.0 - 1.0e-6) { // Near parabolic
		// SGP4 usually caps this, but if we wanted to error:
		// return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonEccentricityTooHigh, Value: e}
		e = 1.0 - 1.0e-6
	}

	// Call CalculateFinalPositionVelocity equivalent
	// This is where short-period perturbations are applied.
	// Constants cosio, sinio, x3thm1, x1mth2, xlcof, aycof are from the *initial* inclination (elems.cosio etc)
	// x7thm1 was calculated during Initialize as well
	beta2 := 1.0 - e*e // current beta_sq, based on current 'e'
	if beta2 < 0.0 {   // Should not happen if e is capped correctly to < 1
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonBeta2Negative, Value: beta2, Message: fmt.Sprintf("based on current eccentricity e=%.6e", e)}
	}
	xn := xke / math.Pow(a, 1.5) // current mean motion (rad/min)

	// Long period periodics (LPP) affecting argument of latitude
	axn := e * math.Cos(omega) // LPP term for L' (e_k * cos(omega_k))
	temp11_lpp := 1.0 / (a * beta2)
	xll_lpp := temp11_lpp * elems.xlcof * axn
	aynl_lpp := temp11_lpp * elems.aycof
	xlt_lpp := xl + xll_lpp                 // L' = L + Lpp
	ayn_lpp := e*math.Sin(omega) + aynl_lpp // LPP term for L' (e_k * sin(omega_k))

	elsq := axn*axn + ayn_lpp*ayn_lpp // (e_k)^2, where e_k is the LPP-perturbed eccentricity vector magnitude
	if elsq >= 1.0 {
		// This can happen if perturbations drive eccentricity too high.
		// For robust code, might cap elsq or return error.
		// libsgp4 throws "Error: (elsq >= 1.0)"
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonPerturbedEccSqTooHigh, Value: elsq}
	}

	// Solve Kepler's equation for L' (eccentric longitude)
	// capu is M' = L' - Omega - omega (where Omega and omega are current nodek and omegak)
	// In SGP4, this is M_k = L' - Omega_k (argument of latitude from ascending node)
	capu := math.Mod(xlt_lpp-xnode, twoPi) // E_k_mean = L' - Omega_k
	epw := capu                            // Initial guess for eccentric anomaly (perturbed) E_k_ecc

	var sinepw, cosepw, ecose, esine float64
	max_newton_raphson := 1.25 * math.Abs(math.Sqrt(elsq)) // Cap on N-R step

	for i := range 10 {
		sinepw = math.Sin(epw)
		cosepw = math.Cos(epw)
		ecose = axn*cosepw + ayn_lpp*sinepw
		esine = axn*sinepw - ayn_lpp*cosepw

		f_kepler := capu - epw + esine
		if math.Abs(f_kepler) < 1.0e-12 {
			break
		}

		fdot_kepler := 1.0 - ecose
		delta_epw := f_kepler / fdot_kepler

		if i == 0 { // First iteration, apply cap
			if delta_epw > max_newton_raphson {
				delta_epw = max_newton_raphson
			} else if delta_epw < -max_newton_raphson {
				delta_epw = -max_newton_raphson
			}
		} else {
			// Second-order Newton-Raphson correction (matches libsgp4)
			delta_epw = f_kepler / (fdot_kepler + 0.5*esine*delta_epw)
		}
		epw += delta_epw
	}

	// Short period preliminary quantities
	temp21_sp := max(1.0-elsq, 0.0) // Ensure non-negative for sqrt
	pl := a * temp21_sp             // semi-latus rectum p_k = a_k * (1 - (e_k)^2)
	if pl < 0.0 {                   // Should be caught by elsq >= 1.0 or a < 0 (if tempa becomes very negative)
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonSemiLatusRectumNegative, Value: pl}
	}

	r_val := a * (1.0 - ecose) // distance from primary focus r_k = a_k * (1 - e_k cos(E_k_ecc))
	if r_val == 0.0 {
		r_val = 1e-9
	} // Avoid division by zero
	temp31_sp := 1.0 / r_val
	rdot_val := xke * math.Sqrt(a) * esine * temp31_sp // r_dot_k
	rfdot_val := xke * math.Sqrt(pl) * temp31_sp       // r_k * f_dot_k (f is true anomaly)

	temp32_sp := a * temp31_sp       // a_k / r_k
	betal_sp := math.Sqrt(temp21_sp) // beta_k = sqrt(1 - (e_k)^2)
	temp33_sp := 0.0
	if (1.0 + betal_sp) != 0.0 {
		temp33_sp = 1.0 / (1.0 + betal_sp)
	} else {
		// This case (betal_sp = -1) should not happen for real eccentricities
		temp33_sp = 1.0e12 // Avoid division by zero, effectively making next terms small
	}

	cosu_sp := temp32_sp * (cosepw - axn + ayn_lpp*esine*temp33_sp)
	sinu_sp := temp32_sp * (sinepw - ayn_lpp - axn*esine*temp33_sp)
	u_sp := math.Atan2(sinu_sp, cosu_sp)

	sin2u_sp := 2.0 * sinu_sp * cosu_sp
	cos2u_sp := 2.0*cosu_sp*cosu_sp - 1.0

	// Short period perturbations (SPP)
	// Constants x3thm1, x1mth2, cosio, sinio, x7thm1 are from initial elements
	temp41_spp := 0.0
	if pl != 0.0 {
		temp41_spp = 1.0 / pl
	} else {
		temp41_spp = 1.0e12
	}

	temp42_spp := ck2 * temp41_spp
	temp43_spp := temp42_spp * temp41_spp

	rk_spp := r_val*(1.0-1.5*temp43_spp*betal_sp*elems.x3thm1) + 0.5*temp42_spp*elems.x1mth2*cos2u_sp
	uk_spp := u_sp - 0.25*temp43_spp*elems.x7thm1*sin2u_sp
	xnodek_spp := xnode + 1.5*temp43_spp*elems.cosio*sin2u_sp
	xinck_spp := elems.incl + 1.5*temp43_spp*elems.cosio*elems.sinio*cos2u_sp
	rdotk_spp := rdot_val - xn*temp42_spp*elems.x1mth2*sin2u_sp
	rfdotk_spp := rfdot_val + xn*temp42_spp*(elems.x1mth2*cos2u_sp+1.5*elems.x3thm1)

	// Orientation vectors (using perturbed elements)
	sinuk := math.Sin(uk_spp)
	cosuk := math.Cos(uk_spp)
	sinik := math.Sin(xinck_spp)
	cosik := math.Cos(xinck_spp)
	sinnok := math.Sin(xnodek_spp)
	cosnok := math.Cos(xnodek_spp)

	xmx := -sinnok * cosik
	xmy := cosnok * cosik

	// Position in ECI (km)
	// ux, uy, uz are components of position unit vector in ECI
	ux := xmx*sinuk + cosnok*cosuk
	uy := xmy*sinuk + sinnok*cosuk
	uz := sinik * sinuk
	posX := rk_spp * ux * xkmper
	posY := rk_spp * uy * xkmper
	posZ := rk_spp * uz * xkmper

	// Velocity in ECI (km/s)
	// vx, vy, vz are components of (velocity / (r_k * f_dot_k)) unit vector, or similar
	vx_orient := xmx*cosuk - cosnok*sinuk
	vy_orient := xmy*cosuk - sinnok*sinuk
	vz_orient := sinik * cosuk

	// rdotk_spp is in ER/min, rfdotk_spp is in ER/min
	// Convert to km/s by multiplying by (xkmper / 60.0)
	vFactor := xkmper / 60.0
	velX := (rdotk_spp*ux + rfdotk_spp*vx_orient) * vFactor
	velY := (rdotk_spp*uy + rfdotk_spp*vy_orient) * vFactor
	velZ := (rdotk_spp*uz + rfdotk_spp*vz_orient) * vFactor

	// Check for decay (physical condition)
	if rk_spp < 1.0 {
		// Satellite has decayed. SGP4 docs say prediction is unreliable.
		return Eci{}, &SatelliteDecayedError{Tsince: tsince, Radius: rk_spp}
	}

	return Eci{
		DateTime: tle.EpochTime().Add(time.Duration(tsince * float64(time.Minute))),
		Position: Vector{X: posX, Y: posY, Z: posZ},
		Velocity: Vector{X: velX, Y: velY, Z: velZ},
	}, nil
}

// GreenwichSiderealTime calculates the Greenwich Mean Sidereal Time
func (eci *Eci) GreenwichSiderealTime() float64 {
	jd := julianDateTime(eci.DateTime)
	t := (jd - 2451545.0) / 36525.0

	gmst_deg := 280.46061837 +
		360.98564736629*(jd-2451545.0) +
		0.000387933*t*t -
		t*t*t/38710000.0

	gmst_deg = math.Mod(gmst_deg, 360.0)
	if gmst_deg < 0 {
		gmst_deg += 360.0
	}
	return gmst_deg * deg2rad
}

// FindPositionAtTime propagates the TLE to a specific absolute time.
func (tle *TLE) FindPositionAtTime(t time.Time) (Eci, error) {
	// Calculate tsince (time since TLE epoch in minutes)
	tsince := t.Sub(tle.EpochTime()).Minutes()
	return tle.FindPosition(tsince)
}

// GeneratePasses predicts satellite passes over a ground station within a given time window.
// lat, lng are observer's geodetic latitude/longitude in degrees.
// alt is observer's altitude in meters above sea level.
// start, stop are the time window boundaries.
// step is the time step for propagation.
func (tle *TLE) GeneratePasses(obsLat, obsLng, obsAltMeters float64, start, stop time.Time, step time.Duration) ([]PassDetails, error) {
	return tle.GeneratePassesWithHorizon(&Horizon{}, obsLat, obsLng, obsAltMeters, start, stop, step)
}

// GeneratePasses predicts satellite passes over a ground station within a given time window and with a specified horizon profile.
// horizon defines the elevation angle threshold for AOS and LOS, allowing for custom horizon profiles (e.g., obstructions).
// lat, lng are observer's geodetic latitude/longitude in degrees.
// alt is observer's altitude in meters above sea level.
// start, stop are the time window boundaries.
// step is the time step for propagation.
func (tle *TLE) GeneratePassesWithHorizon(horizon *Horizon, obsLat, obsLng, obsAltMeters float64, start, stop time.Time, step time.Duration) ([]PassDetails, error) {
	if start.After(stop) {
		return nil, fmt.Errorf("start time must be before stop time")
	}

	if step <= 0 {
		return nil, fmt.Errorf("step duration must be positive")
	}

	observer := &Location{
		Latitude:  obsLat,
		Longitude: obsLng,
		Altitude:  obsAltMeters,
	}

	var passes []PassDetails

	getLookAngleAtTime := func(t time.Time) (*Observation, error) {
		eciState, err := tle.FindPositionAtTime(t)
		if err != nil {
			return nil, err
		}
		sv := &StateVector{
			X: eciState.Position.X, Y: eciState.Position.Y, Z: eciState.Position.Z,
			VX: eciState.Velocity.X, VY: eciState.Velocity.Y, VZ: eciState.Velocity.Z,
		}
		return sv.GetLookAngle(observer, t)
	}

	// findCrossingPoint finds the exact second the satellite crosses the horizon (0° elevation)
	// between time1 and time2. `findingAOS` should be true for AOS, false for LOS.
	findCrossingPoint := func(horizon *Horizon, time1, time2 time.Time, findingAOS bool) (time.Time, error) {
		var middleTime time.Time
		var err error
		var newObserverAngle *Observation

		/*
		 * loop until we zeroed in on the root below a second time diff
		 */
		for time2.Sub(time1).Seconds() > 1.0 {
			middleTime = time1.Add(time.Second * time.Duration(time2.Sub(time1).Seconds()/2.0))

			/*
			 * calculate elevation at time
			 */
			newObserverAngle, err = getLookAngleAtTime(middleTime)
			if err != nil {
				return time.Time{}, err
			}
			if (newObserverAngle.LookAngles.Elevation > horizon.GetElevation(newObserverAngle.LookAngles.Azimuth)) == findingAOS {
				time2 = middleTime
			} else {
				time1 = middleTime
			}
		}

		/*
		 * go back/forward 1second until below the horizon
		 */
		for newObserverAngle.LookAngles.Elevation > horizon.GetElevation(newObserverAngle.LookAngles.Azimuth) {
			middleTime = middleTime.Add(time.Second * time.Duration(func() int {
				if findingAOS {
					return -1
				}
				return 1
			}()))
			newObserverAngle, err = getLookAngleAtTime(middleTime)
			if err != nil {
				return time.Time{}, err
			}
		}

		return middleTime, nil
	}

	// findMaxElevation finds the maximum elevation between aos and los times.
	findMaxElevation := func(aos, los time.Time) (maxEl float64, maxElTime time.Time, err error) {
		timeStep := (los.Sub(aos).Seconds()) / 9.0
		currentTime := aos       //! current time
		endTime := los           //! end time of search period
		var maxElevation float64 //! max elevation
		var maxElevationTime time.Time

		for timeStep > 1.0 {
			maxElevation = -math.MaxFloat64

			for currentTime.Before(endTime) {
				/*
				 * calculate elevation at time
				 */
				angle, err := getLookAngleAtTime(currentTime)
				if err != nil {
					return math.NaN(), time.Time{}, err
				}
				if angle.LookAngles.Elevation > maxElevation {
					/*
					 * still going up
					 */
					maxElevation = angle.LookAngles.Elevation
					maxElevationTime = currentTime
					/*
					 * move time along
					 */
					currentTime = currentTime.Add(time.Second * time.Duration(timeStep))
					if currentTime.After(endTime) {
						/*
						 * dont go past end time
						 */
						currentTime = endTime
					}
				} else {
					/*
					 * stop
					 */
					break
				}
			}
			/*
			 * make end time to current time
			 */
			endTime = currentTime
			/*
			 * make current time of search period to 2 time steps back
			 */
			currentTime = currentTime.Add(time.Second * time.Duration(-2.0*timeStep))
			/*
			 * recalculate time step
			 */
			timeStep = (endTime.Sub(currentTime).Seconds()) / 9.0
		}
		return maxElevation, maxElevationTime, nil
	}

	// Main pass generation loop
	isGeostationary := tle.IsGeostationary()
	currentTime := start.Round(time.Second) // Round to nearest second for cleaner pass times
	currentObservation, err := getLookAngleAtTime(currentTime)
	if err != nil {
		return nil, fmt.Errorf("error getting initial elevation: %w", err)
	}
	foundAos := false
	foundLos := false
	var aosTime time.Time

	if !isGeostationary {
		for currentObservation.LookAngles.Elevation >= horizon.GetElevation(currentObservation.LookAngles.Azimuth) {
			currentObservation, err = getLookAngleAtTime(currentTime)
			if err != nil {
				break // Stop if we encounter an error during propagation
			}
			currentTime = currentTime.Add(step * -1) // Step backwards in time
		}
	} else if currentObservation.LookAngles.Elevation > horizon.GetElevation(currentObservation.LookAngles.Azimuth) {
		foundAos = true
		aosTime = currentTime
	}

	var passDataPoints []PassDataPoint
	for currentTime.Before(stop) || foundAos {
		var losTime time.Time
		nextTime := currentTime.Add(step)
		/*
		 * calculate next elevation
		 */
		nextObservation, err := getLookAngleAtTime(nextTime)
		nextObservationElevationThreshold := horizon.GetElevation(nextObservation.LookAngles.Azimuth)
		currentObservationElevationThreshold := horizon.GetElevation(currentObservation.LookAngles.Azimuth)
		if err != nil {
			return passes, err
		}
		if currentObservation.LookAngles.Elevation < currentObservationElevationThreshold && nextObservation.LookAngles.Elevation >= nextObservationElevationThreshold {
			/*
			 * find the point at which the satellite crossed the horizon
			 */
			aosTime, err = findCrossingPoint(horizon, currentTime, nextTime, true)
			if err == nil {
				foundAos = true
			}
		} else if isGeostationary && !stop.Before(nextTime) && foundAos {
			/*
			 * For geostationary satellites, if we are still within the time window and have found an AOS,
			 * we can assume LOS is at the end of the time window since they won't set.
			 */
			losTime = stop
			foundLos = true
		} else if currentObservation.LookAngles.Elevation > currentObservationElevationThreshold && nextObservation.LookAngles.Elevation <= nextObservationElevationThreshold && foundAos {
			/*
			 * already have the aos, but now the satellite is below the horizon,
			 * so find the los
			 */
			losTime, err = findCrossingPoint(horizon, currentTime, nextTime, false)
			if err == nil {
				foundLos = true
			} else {
				foundAos = false // Reset AOS if we fail to find LOS, to avoid creating a pass with only AOS
			}
		}
		if foundAos && foundLos {
			maxElevation, maxElevationTime, _ := findMaxElevation(aosTime, losTime)
			aosObs, _ := getLookAngleAtTime(aosTime)
			maxElObs, _ := getLookAngleAtTime(maxElevationTime)
			losObs, _ := getLookAngleAtTime(losTime)

			passDataPoints = append(passDataPoints,
				passDataPointFromObservation(aosTime, aosObs),
				passDataPointFromObservation(maxElevationTime, maxElObs),
				passDataPointFromObservation(losTime, losObs),
			)
			slices.SortFunc(passDataPoints, func(a, b PassDataPoint) int {
				return int(a.Timestamp.Sub(b.Timestamp).Seconds())
			})

			pd := PassDetails{
				AOS:              aosTime,
				LOS:              losTime,
				AOSAzimuth:       aosObs.LookAngles.Azimuth,
				LOSAzimuth:       losObs.LookAngles.Azimuth,
				MaxElevation:     maxElevation,
				MaxElevationAz:   maxElObs.LookAngles.Azimuth,
				MaxElevationTime: maxElevationTime,
				AOSObservation:   *aosObs,
				LOSObservation:   *losObs,
				MaxElObservation: *maxElObs,
				Duration:         losTime.Sub(aosTime),
				DataPoints:       passDataPoints,
			}

			passes = append(passes, pd)
			foundAos = false
			foundLos = false
			passDataPoints = []PassDataPoint{}
			nextTime = losTime.Add(time.Duration(30) * time.Minute) // Skip 30 minutes
		}

		currentTime = nextTime
		currentObservation = nextObservation

		if foundAos && !foundLos {
			passDataPoints = append(passDataPoints, passDataPointFromObservation(currentTime, currentObservation))
		}
	}

	return passes, nil
}

func passDataPointFromObservation(time time.Time, obs *Observation) PassDataPoint {
	return PassDataPoint{
		Timestamp: time,
		Azimuth:   obs.LookAngles.Azimuth,
		Elevation: obs.LookAngles.Elevation,
		Range:     obs.LookAngles.Range,
		RangeRate: obs.LookAngles.RangeRate,
	}
}

func julianDateTime(t_utc time.Time) float64 {
	y := float64(t_utc.Year())
	m := float64(t_utc.Month())
	d := float64(t_utc.Day())
	h := float64(t_utc.Hour())
	min := float64(t_utc.Minute())
	s := float64(t_utc.Second())
	ns := float64(t_utc.Nanosecond())

	if m <= 2 {
		y--
		m += 12
	}
	A_jd := math.Floor(y / 100.0)
	B_jd := 2 - A_jd + math.Floor(A_jd/4.0)
	JD_day := math.Floor(365.25*(y+4716.0)) +
		math.Floor(30.6001*(m+1.0)) +
		d + B_jd - 1524.5
	dayFrac := (h + min/60.0 + s/3600.0 + ns/3.6e12) / 24.0
	return JD_day + dayFrac
}

func wrapLongitude(lon float64) float64 {
	lon = math.Mod(lon, twoPi)
	if lon > math.Pi {
		lon -= twoPi
	} else if lon < -math.Pi {
		lon += twoPi
	}
	return lon
}

// AcTan calculates the arctangent in the correct quadrant based on sin and cos values.
// It replicates the behavior of the C++ libsgp4 Util::AcTan function.
func AcTan(sinx, cosx float64) float64 {
	// Handle the case where cosx is effectively zero to avoid division by zero
	// and determine the correct quadrant based on the sign of sinx.
	if math.Abs(cosx) < 1e-14 { // Using a small epsilon, similar to other checks
		if sinx > 0.0 {
			return math.Pi / 2.0
		} else {
			// Note: C++ returns 3*PI/2. This is equivalent to -PI/2 in terms of angle,
			// but let's stick exactly to the C++ logic for maximum fidelity.
			// If the downstream code expects [-PI, PI] range, you might consider
			// returning -math.Pi / 2.0 instead, but C++ uses [0, 2*PI).
			return 3.0 * math.Pi / 2.0
		}
	} else {
		// If cosx is positive, the angle is in quadrants I or IV.
		// math.Atan will return the correct angle in [-PI/2, PI/2].
		if cosx > 0.0 {
			return math.Atan(sinx / cosx)
		} else {
			// If cosx is negative, the angle is in quadrants II or III.
			// math.Atan will return an angle in [-PI/2, PI/2], but the actual
			// angle is in [PI/2, 3*PI/2]. Adding PI corrects the quadrant.
			// This also correctly handles the case where sinx is zero and cosx is negative
			// (angle is PI).
			return math.Pi + math.Atan(sinx/cosx) // Adding PI shifts to correct quadrant
		}
	}
}
