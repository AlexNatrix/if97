package secondRegionMeta

import (
	"math"

	"if97.com/cmd/lib/utils/constants"
)

var (
	Jno = [][]float64{
		{0, -00.96937268393049e1},
		{1, 000.10087275970006e2},
		{-5, -0.56087911283020e-2},
		{-4, 00.71452738081455e-1},
		{-3, -0.40710498223928},
		{-2, 00.14240819171444e1},
		{-1, -0.43839511319450e1},
		{2, -00.28408632460772},
		{3, 000.21268463753307e-1}}
	IJnr = [][]float64{
		{1, 0, -0.73362260186506e-2},
		{1, 2, -0.88223831943146e-1},
		{1, 5, -0.72334555213245e-1},
		{1, 11, -.40813178534455e-2},
		{2, 1, 00.20097803380207e-2},
		{2, 7, -0.53045921898642e-1},
		{2, 16, -.76190409086970e-2},
		{3, 4, -0.63498037657313e-2},
		{3, 16, -.86043093028588e-1},
		{4, 7, 00.75321581522770e-2},
		{4, 10, -.79238375446139e-2},
		{5, 9, -0.22888160778447e-3},
		{5, 10, -.26456501482810e-2}}
)

const (
	Tref = 540
	pRef = 1
	Name = "Region 2 metastable-vapour"
)


type Region struct {
	Name string
}

var REGION2META = Region{
	Name: Name,
}

func (r *Region)GetName()string{
	return r.Name
}


/**
 * Dimensionless Gibbs free energy.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return
 */
func gammaO(pi float64, tau float64) float64 {

	var out float64 = math.Log(pi)

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

	return 1 / pi
}

/**
 * Second partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @return
 */
func gammaOPiPi(pi float64) float64 {

	return -1 / (pi * pi)
}

/**
 * Second partial derivative with respect to pi & tau.
 *
 * @return
 */
func gammaOPiTau(...float64) float64 {

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
		out += ijnr[2] * math.Pow(pi, ijnr[0]) * math.Pow(tau-0.5, ijnr[1])
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
		out += ijnr[2] * ijnr[0] * math.Pow(pi, ijnr[0]-1) * math.Pow(tau-0.5, ijnr[1])
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
		out += ijnr[2] * ijnr[0] * (ijnr[0] - 1) * math.Pow(pi, ijnr[0]-2) * math.Pow(tau-0.5, ijnr[1])
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
		out += ijnr[2] * ijnr[0] * math.Pow(pi, ijnr[0]-1) * ijnr[1] * math.Pow(tau-0.5, ijnr[1]-1)
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
		out += ijnr[2] * math.Pow(pi, ijnr[0]) * ijnr[1] * math.Pow(tau-0.5, ijnr[1]-1)
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
		out += ijnr[2] * math.Pow(pi, ijnr[0]) * ijnr[1] * (ijnr[1] - 1) * math.Pow(tau-0.5, ijnr[1]-2)
	}
	return out
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
		(1 - pi*pi*gammaRPiPi(pi, tau) + x*x/
			(tau*tau*(gammaOTauTau(tau)+gammaRTauTau(pi, tau)))) * 1e3 * constants.R * temperature)
}
