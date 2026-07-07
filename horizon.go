package sgp4

import (
	"math"
)

type Horizon [360]float64

func (h *Horizon) GetElevation(azimuth float64) float64 {
	for azimuth < 0 {
		azimuth += 360
	}
	azimuth = math.Floor(azimuth)
	azimuth = math.Mod(azimuth, 360)
	return h[int(azimuth)]
}
