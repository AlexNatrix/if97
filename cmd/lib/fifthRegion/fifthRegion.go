package fifthRegion

import (
	"math"

	"if97.com/cmd/lib/utils/constants"
)

const (
	Name = "Region 5"
	Tref = 1000
	pRef = 1
)

var (
	Jno = [][]float64{
		{00, -.131799836742010e2},
		{01, 0.685408416344340e1},
		{-3, -.248051489334660e-1},
		{-2, 0.369015349803330},
		{-1, -.311613182139250e1},
		{02, -.329616265389170}}
	IJnr = [][]float64{
		{1, 1, 0.15736404855259e-2},
		{1, 2, 0.90153761673944e-3},
		{1, 3, -.50270077677648e-2},
		{2, 3, 0.22440037409485e-5},
		{2, 9, -.41163275453471e-5},
		{3, 7, 0.37919454822955e-7}}
)

type Region struct {
	Name string
}

var REGION5 = Region{
	Name: Name,
}

func (r *Region)GetName() string{
	return r.Name
}

type SubRegion int64

const (
	_ SubRegion = iota
	A
	B
	C
)

func enthalpy2bc(pressure float64) float64 {

	n := []float64{
		0.90584278514712e3,
		-.67955786399241,
		0.12809002730136e-3,
		0.26526571908428e4,
		0.45257578905948e1}

	return n[3] + math.Sqrt((pressure-n[4])/n[2])
}

/**
 * Dimensionless Gibbs free energy.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaO(pi float64, tau float64) float64 {

	out := math.Log(pi)

	for _, jno := range Jno {
		out += jno[1] * math.Pow(tau, jno[0])
	}
	return out
}

/**
 * First partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @return
 */
func gammaOPi(pi float64) float64 {
	return 1.0 / pi
}

/**
 * Second partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @return
 */
func gammaOPiPi(pi float64) float64 {
	return -1.0 / (pi * pi)
}

/**
 * Second partial derivative with respect to pi & tau.
 *
 * @return
 */
func gammaOPiTau() float64 {
	return 0
}

/**
 * First partial derivative with respect to tau.
 *
 * @param tau dimensionless temperature
 * @return
 */
func gammaOTau(tau float64) float64 {

	var out float64 = 0

	for _, jno := range Jno {
		out += jno[1] * jno[0] * math.Pow(tau, jno[0]-1)
	}
	return out
}

/**
 * Second partial derivative with respect to tau.
 *
 * @param tau dimensionless temperature
 * @return
 */
func gammaOTauTau(tau float64) float64 {

	var out float64 = 0

	for _, jno := range Jno {
		out += jno[1] * jno[0] * (jno[0] - 1) * math.Pow(tau, jno[0]-2)
	}
	return out
}

/**
 * Dimensionless Gibbs free energy.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaR(pi float64, tau float64) float64 {

	var out float64 = 0

	for _, ijnr := range IJnr {
		out += ijnr[2] * math.Pow(pi, ijnr[0]) * math.Pow(tau, ijnr[1])
	}
	return out
}

/**
 * First partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaRPi(pi float64, tau float64) float64 {

	var out float64 = 0

	for _, ijnr := range IJnr {
		out += ijnr[2] * ijnr[0] * math.Pow(pi, ijnr[0]-1) * math.Pow(tau, ijnr[1])
	}
	return out
}

/**
 * Second partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaRPiPi(pi float64, tau float64) float64 {

	var out float64 = 0

	for _, ijnr := range IJnr {
		out += ijnr[2] * ijnr[0] * (ijnr[0] - 1) * math.Pow(pi, ijnr[0]-2) * math.Pow(tau, ijnr[1])
	}
	return out
}

/**
 * Second partial derivative with respect to pi & tau.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaRPiTau(pi float64, tau float64) float64 {

	var out float64 = 0

	for _, ijnr := range IJnr {
		out += ijnr[2] * ijnr[0] * math.Pow(pi, ijnr[0]-1) * ijnr[1] * math.Pow(tau, ijnr[1]-1)
	}
	return out
}

/**
 * First partial derivative with respect to tau.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaRTau(pi float64, tau float64) float64 {

	var out float64 = 0

	for _, ijnr := range IJnr {
		out += ijnr[2] * math.Pow(pi, ijnr[0]) * ijnr[1] * math.Pow(tau, ijnr[1]-1)
	}
	return out
}

/**
 * Second partial derivative with respect to tau.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaRTauTau(pi float64, tau float64) float64 {

	var out float64 = 0

	for _, ijnr := range IJnr {
		out += ijnr[2] * math.Pow(pi, ijnr[0]) * ijnr[1] * (ijnr[1] - 1) * math.Pow(tau, ijnr[1]-2)
	}
	return out
}

func GetSubRegion(pressure float64, enthalpy float64) SubRegion {
	if pressure > 4 {
		if enthalpy < enthalpy2bc(pressure) {
			return C
		} else {
			return B
		}

	} else {
		return A
	}

}


func (r *Region)HeatCapacityRatioPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	x := 1 + pi*gammaRPi(pi, tau) - tau*pi*gammaRPiTau(pi, tau)

	return 1 / (1 + x*x/(tau*tau*(gammaOTauTau(tau)+gammaRTauTau(pi, tau))*(1-pi*pi*gammaRPiPi(pi, tau))))
}

func  (r *Region)SpecificVolumePH(p float64, h float64) float64 {

	T := r.TemperaturePH(p, h)

	return r.SpecificVolumePT(p, T)
}


func (r *Region)SpecificVolumePS(p float64, s float64) float64 {

	T := r.TemperaturePS(p, s)

	return r.SpecificVolumePT(p, T)
}

func (r *Region)IsentropicExponentPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	gRPi := gammaRPi(pi, tau)
	x := 1 + pi*gRPi
	y := x * tau * tau * (gammaOTauTau(tau) + gammaRTauTau(pi, tau))
	z := x - pi*tau*gammaRPiTau(pi, tau)

	return x * y / ((1-pi*pi*gammaRPiPi(pi, tau))*y + pi*(gammaOPi(pi)+gRPi)*z*z)
}


func (r *Region)IsobaricCubicExpansionCoefficientPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (1 - tau*pi*gammaRPiTau(pi, tau)/(1+pi*gammaRPi(pi, tau))) / temperature
}


func (r *Region)IsothermalCompressibilityPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (1 - pi*pi*gammaRPiPi(pi, tau)) / (1 + pi*gammaRPi(pi, tau)) / pressure
}


func (r *Region)SpecificEnthalpyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return tau * (gammaOTau(tau) + gammaRTau(pi, tau)) * constants.R * temperature
}


func (r *Region)SpecificEntropyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (tau*(gammaOTau(tau)+gammaRTau(pi, tau)) - (gammaO(pi, tau) + gammaR(pi, tau))) * constants.R
}


func (r *Region)SpecificInternalEnergyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (tau*(gammaOTau(tau)+gammaRTau(pi, tau)) - pi*(gammaOPi(pi)+gammaRPi(pi, tau))) * constants.R * temperature
}

func (r *Region)SpecificGibbsFreeEnergyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (gammaO(pi, tau) + gammaR(pi, tau)) * constants.R * temperature
}


func (r *Region)SpecificIsobaricHeatCapacityPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return -tau * tau * (gammaOTauTau(tau) + gammaRTauTau(pi, tau)) * constants.R
}


func (r *Region)SpecificIsochoricHeatCapacityPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	x := 1 + pi*gammaRPi(pi, tau) - tau*pi*gammaRPiTau(pi, tau)

	return (-tau*tau*(gammaOTauTau(tau)+gammaRTauTau(pi, tau)) - x*x/
		(1-pi*pi*gammaRPiPi(pi, tau))) * constants.R
}


func (r *Region)SpecificVolumePT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return pi * (gammaOPi(pi) + gammaRPi(pi, tau)) / 1e3 * constants.R * temperature / pressure
}


func (r *Region)SpeedOfSoundPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	gRPi := gammaRPi(pi, tau)
	x := 1 + pi*gRPi - tau*pi*gammaRPiTau(pi, tau)

	return math.Sqrt((1 + 2*pi*gRPi + pi*pi*gRPi*gRPi) /
		(1 - pi*pi*gammaRPiPi(pi, tau) + x*x/(tau*tau*(gammaOTauTau(tau)+gammaRTauTau(pi, tau)))) * 1e3 * constants.R * temperature)
}


func (r *Region)TemperaturePH(pressure float64, enthalpy float64) float64 {
	return math.NaN()
}


func (r *Region)PressureHS(h float64, s float64) float64 {
	return math.NaN()
}


func (r *Region)SpecificEntropyRhoT(rho float64, T float64) float64 {
	return math.NaN()
}


func (r *Region)TemperatureHS(h float64, s float64) float64 {
	return math.NaN()
}


func (r *Region)TemperaturePS(p float64, s float64) float64 {
	return math.NaN()
}
