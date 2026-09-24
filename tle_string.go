package sgp4

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// TLEString formats a TLE as two fixed-width lines, preceded by the name if set.
// The name is limited to 24 characters. Checksums are calculated from the
// formatted lines rather than the stored fields.
func (tle *TLE) TLEString() (string, error) {
	if tle == nil {
		return "", fmt.Errorf("cannot format a nil TLE")
	}

	id, err := tle.NoradIDString()
	if err != nil {
		return "", err
	}

	dot, err := formatTLEDot(tle.MeanMotionDot)
	if err != nil {
		return "", fmt.Errorf("mean motion dot: %w", err)
	}
	dot2 := formatTLEExponent(tle.MeanMotionDot2)
	bstar := formatTLEExponent(tle.Bstar)

	// TLE stores seven fractional digits without rounding
	eccentricityDigits := strings.TrimPrefix(strconv.FormatFloat(tle.Eccentricity, 'f', -1, 64), "0.")
	eccentricity := (eccentricityDigits + "0000000")[:7]

	line1 := fmt.Sprintf("1 %s%c %-8s %02d%012.8f %s %s %s 0 %4d",
		id, tle.Classification, tle.International, tle.EpochYear%100, tle.EpochDay,
		dot, dot2, bstar, tle.ElementNumber)
	line2 := fmt.Sprintf("2 %s %8.4f %8.4f %s %8.4f %8.4f %11.8f%5d",
		id, tle.Inclination, tle.RightAscension, eccentricity, tle.ArgOfPerigee, tle.MeanAnomaly, tle.MeanMotion,
		tle.RevolutionNumber)
	if len(line1) != 68 || len(line2) != 68 {
		return "", fmt.Errorf("formatted TLE lines have invalid lengths: %d and %d", len(line1), len(line2))
	}
	check1, _ := CalculateChecksum(line1 + "0")
	check2, _ := CalculateChecksum(line2 + "0")
	result := fmt.Sprintf("%s%d\n%s%d", line1, check1, line2, check2)
	name := strings.TrimSpace(tle.Name)
	if name != "" {
		if len(name) > 24 {
			name = string([]rune(name)[:24])
		}
		result = string(name) + "\n" + result
	}
	return result, nil
}

func formatTLEDot(value float64) (string, error) {
	formatted := fmt.Sprintf("% .8f", value)
	if len(formatted) != 11 || formatted[1] != '0' {
		return "", fmt.Errorf("value cannot fit the TLE field: %v", value)
	}
	return formatted[:1] + formatted[2:], nil
}

func formatTLEExponent(value float64) string {
	exponent := 0
	if value != 0 {
		exponent = 1 + int(math.Floor(math.Log10(math.Abs(value)))) // The TLE mantissa is 0.xxxxx, not 1.xxxx.
	}
	mant := int(math.Round(value / math.Pow(10, float64(exponent-5))))
	return fmt.Sprintf("% 06d%+2d", mant, exponent)
}
