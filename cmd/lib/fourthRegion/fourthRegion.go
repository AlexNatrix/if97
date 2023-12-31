package fourthRegion

import (
	"math"

	"if97.com/cmd/lib/firstRegion"
	rangeError "if97.com/cmd/lib/region/range_error"
	"if97.com/cmd/lib/secondRegion"
	"if97.com/cmd/lib/thirdRegion"

	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
)

const (
	Name            = "Region 4"
	Tref            = 1
	pRef            = 1
	ITERATION_LIMIT = 100
	TOLERANCE       = 1e-9
)

var (
	n = []float64{
		00.11670521452767e4,
		-0.72421316703206e6,
		-0.17073846940092e2,
		00.12020824702470e5,
		-0.32325550322333e7,
		00.14915108613530e2,
		-0.48232657361591e4,
		00.40511340542057e6,
		-0.23855557567849,
		00.65017534844798e3}
	IJnH = [][]float64{
		{0, 0, 0.600073641753024},
		{1, 1, -.936203654849857e1},
		{1, 3, 0.246590798594147e2},
		{1, 4, -.107014222858224e3},
		{1, 36, -.915821315805768e14},
		{5, 03, -.862332011700662e4},
		{7, 00, -.235837344740032e2},
		{8, 24, 0.252304969384128e18},
		{14, 16, -.389718771997719e19},
		{20, 16, -.333775713645296e23},
		{22, 03, 0.356499469636328e11},
		{24, 18, -.148547544720641e27},
		{28, 8, 0.330611514838798e19},
		{36, 24, 0.813641294467829e38}}
	IJnS = [][]float64{
		{0, 0, .639767553612785},
		{1, 1, -.129727445396014e2},
		{1, 32, -.224595125848403e16},
		{4, 7, .177466741801846e7},
		{12, 4, .717079349571538e10},
		{12, 14, -.378829107169011e18},
		{16, 36, -.955586736431328e35},
		{24, 10, .187269814676188e24},
		{28, 0, .119254746466473e12},
		{32, 18, .110649277244882e37}}
	IJnHS = [][]float64{
		{0, 0, .179882673606601},
		{0, 3, -.267507455199603},
		{0, 12, .116276722612600e1},
		{1, 0, .147545428713616},
		{1, 1, -.512871635973248},
		{1, 2, .421333567697984},
		{1, 5, .563749522189870},
		{2, 0, .429274443819153},
		{2, 5, -.335704552142140e1},
		{2, 8, .108890916499278e2},
		{3, 0, -.248483390456012},
		{3, 2, .304153221906390},
		{3, 3, -.494819763939905},
		{3, 4, .107551674933261e1},
		{4, 0, .733888415457688e-1},
		{4, 1, .140170545411085e-1},
		{5, 1, -.106110975998808},
		{5, 2, .168324361811875e-1},
		{5, 4, .125028363714877e1},
		{5, 16, .101316840309509e4},
		{6, 6, -.151791558000712e1},
		{6, 8, .524277865990866e2},
		{6, 22, .230495545563912e5},
		{8, 1, .249459806365456e-1},
		{10, 20, .210796467412137e7},
		{10, 36, .366836848613065e9},
		{12, 24, -.144814105365163e9},
		{14, 1, -.179276373003590e-2},
		{14, 28, .489955602100459e10},
		{16, 12, .471262212070518e3},
		{16, 32, -.829294390198652e11},
		{18, 14, -.171545662263191e4},
		{18, 22, .355777682973575e7},
		{18, 36, .586062760258436e12},
		{20, 24, -.129887635078195e8},
		{28, 36, .317247449371057e11}}
)

var (
	p0   = REGION4.SaturationPressureT(constants.T0)
	s2   = secondRegion.REGION2.SpecificEntropyPT(p0, constants.T0)
	ps13 = REGION4.SaturationPressureT(constants.T13) // (16.529 MPa) [MPa]
	hs13 = firstRegion.REGION1.SpecificEnthalpyPT(ps13, constants.T13)
	ss13 = firstRegion.REGION1.SpecificEntropyPT(ps13, constants.T13)
	hs23 = secondRegion.REGION2.SpecificEnthalpyPT(ps13, constants.T13)
	ss23 = secondRegion.REGION2.SpecificEntropyPT(ps13, constants.T13)
	sc   = thirdRegion.REGION3.SpecificEntropyRhoT(constants.Rhoc, constants.Tc)
)

type Region struct {
	Name string
}

var REGION4 = Region{
	Name: Name,
}

func (r *Region) GetName() string {
	return r.Name
}

func checkB34S(entropy float64) error {

	if entropy < ss13 {
		rangeError.ErrorFromValue(quantity.S, entropy, ss13)

	} else if entropy > ss23 {
		rangeError.ErrorFromValue(quantity.S, entropy, ss23)
	}
	return nil
}

func specificEnthalpy3a(s float64) float64 {

	var eta float64 = 0
	sigma := s / 3.8
	x := []float64{sigma - 1.09, sigma + 0.366e-4, 1700}
	IJn := [][]float64{
		{0, 1, .822673364673336},
		{0, 4, .181977213534479},
		{0, 10, -.112000260313624e-1},
		{0, 16, -.746778287048033e-3},
		{2, 1, -.179046263257381},
		{3, 36, .424220110836657e-1},
		{4, 3, -.341355823438768},
		{4, 16, -.209881740853565e1},
		{5, 20, -.822477343323596e1},
		{5, 36, -.499684082076008e1},
		{6, 4, .191413958471069},
		{7, 2, .581062241093136e-1},
		{7, 28, -.165505498701029e4},
		{7, 32, .158870443421201e4},
		{10, 14, -.850623535172818e2},
		{10, 32, -.317714386511207e5},
		{10, 36, -.945890406632871e5},
		{32, 0, -.139273847088690e-5},
		{32, 6, .631052532240980}}

	for _, ijn := range IJn {
		eta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return eta * x[2]
}

func specificEnthalpy2c3b(s float64) float64 {

	var eta float64 = 0
	sigma := s / 5.9
	x := []float64{sigma - 1.02, sigma - 0.726}
	IJn := [][]float64{
		{0, 0, .104351280732769e1},
		{0, 3, -.227807912708513e1},
		{0, 4, .180535256723202e1},
		{1, 0, .420440834792042},
		{1, 12, -.105721244834660e6},
		{5, 36, .436911607493884e25},
		{6, 12, -.328032702839753e12},
		{7, 16, -.678686760804270e16},
		{8, 2, .743957464645363e4},
		{8, 20, -.356896445355761e20},
		{12, 32, .167590585186801e32},
		{16, 36, -.355028625419105e38},
		{22, 2, .396611982166538e12},
		{22, 32, -.414716268484468e41},
		{24, 7, .359080103867382e19},
		{36, 20, -.116994334851995e41}}

	for _, ijn := range IJn {
		eta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return math.Pow(eta, 4) * 2800
}

func specificEnthalpy2ab(s float64) float64 {

	var eta float64 = 0
	x := []float64{5.21/s - 0.513, s/9.2 - 0.524, 2800}
	IJn := [][]float64{
		{1, 8, -.524581170928788e3},
		{1, 24, -.926947218142218e7},
		{2, 4, -.237385107491666e3},
		{2, 32, .210770155812776e11},
		{4, 1, -.239494562010986e2},
		{4, 2, .221802480294197e3},
		{7, 7, -.510472533393438e7},
		{8, 5, .124981396109147e7},
		{8, 12, .200008436996201e10},
		{10, 1, -.815158509791035e3},
		{12, 0, -.157612685637523e3},
		{12, 7, -.114200422332791e11},
		{18, 10, .662364680776872e16},
		{20, 12, -.227622818296144e19},
		{24, 32, -.171048081348406e32},
		{28, 8, .660788766938091e16},
		{28, 12, .166320055886021e23},
		{28, 20, -.218003784381501e30},
		{28, 22, -.787276140295618e30},
		{28, 24, .151062329700346e32},
		{32, 2, .795732170300541e7},
		{32, 7, .131957647355347e16},
		{32, 12, -.325097068299140e24},
		{32, 14, -.418600611419248e26},
		{32, 24, .297478906557467e35},
		{36, 10, -.953588761745473e20},
		{36, 12, .166957699620939e25},
		{36, 20, -.175407764869978e33},
		{36, 22, .347581490626396e35},
		{36, 28, -.710971318427851e39}}

	for _, ijn := range IJn {
		eta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return math.Exp(eta) * x[2]
}

func CheckHS(enthalpy float64, entropy float64) error {

	if entropy < ss23 {
		return rangeError.ErrorFromValue(quantity.S, entropy, ss23)

	} else if entropy > s2 {
		return rangeError.ErrorFromValue(quantity.S, entropy, s2)

	} else if entropy <= sc {
		h3a := specificEnthalpy3a(entropy)

		if enthalpy > h3a {
			return rangeError.ErrorFromValue(quantity.H, enthalpy, h3a)
		}
	} else if entropy < constants.S2bc {
		h2c3b := specificEnthalpy2c3b(entropy)

		if enthalpy > h2c3b {
			return rangeError.ErrorFromValue(quantity.H, enthalpy, h2c3b)
		}
	} else {
		h2ab := specificEnthalpy2ab(entropy)

		if enthalpy > h2ab {
			return rangeError.ErrorFromValue(quantity.H, enthalpy, h2ab)
		}
	}
	return nil
}

func CheckP(pressure float64) error {

	if pressure < p0 {
		return rangeError.ErrorFromValue(quantity.P, pressure, p0)

	} else if pressure > constants.Pc {
		return rangeError.ErrorFromValue(quantity.P, pressure, constants.Pc)
	}
	return nil
}

func CheckT(temperature float64) error {

	if temperature < constants.T0 {
		return rangeError.ErrorFromValue(quantity.T, temperature, constants.T0)

	} else if temperature > constants.Tc {
		return rangeError.ErrorFromValue(quantity.T, temperature, constants.Tc)
	}
	return nil
}

func densitiesRegion3(pressure float64, enthalpies []float64) []float64 {
	return []float64{
		1.0 / thirdRegion.REGION3.SpecificVolumePH(pressure, enthalpies[0]),
		1.0 / thirdRegion.REGION3.SpecificVolumePH(pressure, enthalpies[1]),
	}
}

/**
 * Gets the derivative, dT/dp, along the saturation line in SI units,
 * according Clausius-Clapeyron.
 *
 * @param pressure [MPa]
 * @return dT/dp [K/Pa]
 */
func DerivativeP(pressure float64) float64 {

	T := REGION4.SaturationTemperatureP(pressure) // [K]
	h := specificEnthalpiesP(pressure)            // [kJ/kg]
	v := specificVolumesP(pressure)               // [m³/kg]

	return T * (v[1] - v[0]) / (h[1] - h[0]) / 1e3
} //TODO Replace by derivative of elementary function

func (r *Region) HeatCapacityRatioPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) IsentropicExponentPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) IsobaricCubicExpansionCoefficientPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.IsobaricCubicExpansionCoefficientPX(pressure, x)
}

func (r *Region) IsobaricCubicExpansionCoefficientPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) IsobaricCubicExpansionCoefficientPX(pressure float64, vapourFraction float64) float64 {

	T := r.SaturationTemperatureP(pressure)
	var alpha []float64

	if pressure > ps13 {
		/*
		   Region 3
		*/
		h := specificEnthalpiesP(pressure)
		rho := densitiesRegion3(pressure, h)
		alpha[0] = thirdRegion.REGION3.IsobaricCubicExpansionCoefficientRhoT(rho[0], T)
		alpha[1] = thirdRegion.REGION3.IsobaricCubicExpansionCoefficientRhoT(rho[1], T)

	} else {
		/*
		   Regions 1 & 2
		*/
		alpha[0] = firstRegion.REGION1.IsobaricCubicExpansionCoefficientPT(pressure, T)
		alpha[1] = secondRegion.REGION2.IsobaricCubicExpansionCoefficientPT(pressure, T)
	}
	return valueX(vapourFraction, alpha)
}

func (r *Region) IsothermalCompressibilityPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.IsothermalCompressibilityPX(pressure, x)
}

func (r *Region) IsothermalCompressibilityPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) IsothermalCompressibilityPX(pressure float64, vapourFraction float64) float64 {

	T := r.SaturationTemperatureP(pressure)
	var kappaT []float64

	if pressure > ps13 {
		/*
		   Region 3
		*/
		h := specificEnthalpiesP(pressure)
		rho := densitiesRegion3(pressure, h)
		kappaT[0] = thirdRegion.REGION3.IsothermalCompressibilityRhoT(rho[0], T)
		kappaT[1] = thirdRegion.REGION3.IsothermalCompressibilityRhoT(rho[1], T)

	} else {
		/*
		   Regions 1 & 2
		*/
		kappaT[0] = firstRegion.REGION1.IsothermalCompressibilityPT(pressure, T)
		kappaT[1] = secondRegion.REGION2.IsothermalCompressibilityPT(pressure, T)
	}
	return valueX(vapourFraction, kappaT)
}

/**
 * [IF97 Supplementary Release S04, June 2014]
 *
 * @param enthalpy specific enthalpy [kJ/kg]
 * @param entropy specific entropy [kJ/kg-K]
 * @return saturation pressure [MPa]
 */

func (r *Region) PressureHS(enthalpy float64, entropy float64) float64 {

	T := r.TemperatureHS(enthalpy, entropy)

	return r.SaturationPressureT(T)
}

/**
 * Boundary saturation pressure for the boundary between regions 3 and 4.
 *
 * @param enthalpy specific enthalpy [kJ/kg]
 * @return saturation pressure [MPa]
 */
func SaturationPressureB34H(enthalpy float64) float64 {

	eta := enthalpy / 2600.0
	var out float64 = 0

	for _, ijn := range IJnH {
		out += ijn[2] * math.Pow(eta-1.02, ijn[0]) * math.Pow(eta-0.608, ijn[1])
	}
	return out * 22
}

/**
 * Boundary saturation pressure for the boundary between regions 3 and 4,
 * only!
 *
 * @param entropy specific entropy [kJ/(kg K)]
 * @return saturation pressure [MPa]
 * @return RangeError
 */
func saturationPressureB34S(entropy float64) (float64, error) {

	err := checkB34S(entropy)
	if err != nil {
		return -1, err
	}

	sigma := entropy / 5.2
	var pi float64 = 0
	x := []float64{sigma - 1.03, sigma - 0.699}

	for _, ijn := range IJnS {
		pi += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return pi * 22, err
}

/**
 * Saturation pressure.
 *
 * Out-of-range exceptions not thrown because this method is also used for
 * finding regions.
 *
 * @param saturationTemperature saturation temperature [K]
 * @return saturation pressure [MPa]
 */
func (r *Region) SaturationPressureT(saturationTemperature float64) float64 {

	Ts_Tref := saturationTemperature / Tref
	theta := Ts_Tref + n[8]/(Ts_Tref-n[9])
	theta2 := theta * theta
	A := theta2 + n[0]*theta + n[1]
	B := n[2]*theta2 + n[3]*theta + n[4]
	C := n[5]*theta2 + n[6]*theta + n[7]

	return math.Pow(2.0*C/(-B+math.Sqrt(B*B-4*A*C)), 4) * pRef
}

/**
 * Saturation temperature.
 *
 * Out-of-range exceptions not thrown because this method is also used for
 * finding regions.
 *
 * @param saturationPressure saturation pressure [MPa]
 * @return saturation temperature [K]
 */
func (r *Region) SaturationTemperatureP(saturationPressure float64) float64 {

	beta := math.Pow(saturationPressure/pRef, 0.25)
	beta2 := beta * beta
	E := beta2 + n[2]*beta + n[5]
	F := n[0]*beta2 + n[3]*beta + n[6]
	G := n[1]*beta2 + n[4]*beta + n[7]
	D := 2 * G / (-F - math.Sqrt(F*F-4*E*G))
	n9plusD := n[9] + D

	return (n9plusD - math.Sqrt(n9plusD*n9plusD-4*(n[8]+n[9]*D))) / 2.0 * Tref
}

/**
 * Method copied from [Numerical Recipes, 2007].
 *
 * @param a
 * @param b
 * @return
 */
func sign(a float64, b float64) float64 {
	if b >= 0 {
		return a
	}
	return -a // simplified using XOR
}

func specificEnthalpiesP(pressure float64) []float64 {
	return []float64{
		REGION4.SpecificEnthalpySaturatedLiquidP(pressure),
		REGION4.SpecificEnthalpySaturatedVapourP(pressure),
	}
}

func (r *Region) SpecificEnthalpyPS(pressure float64, entropy float64) float64 {

	x := r.VapourFractionPS(pressure, entropy)

	return r.SpecificEnthalpyPX(pressure, x)
}

func (r *Region) SpecificEnthalpyPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

/**
 * Specific enthalpy as a function of pressure & vapour fraction.
 *
 * @param pressure pressure [MPa]
 * @param vapourFraction vapour fraction [-]
 * @return specific enthalpy [kJ/kg]
 */
func (r *Region) SpecificEnthalpyPX(pressure float64, vapourFraction float64) float64 {

	h := specificEnthalpiesP(pressure)

	return valueX(vapourFraction, h)
}

func signum(a float64) float64 {
	switch {
	case a < 0:
		return -1
	case a > 0:
		return +1
	}
	return 0
}

/**
 * Specific enthalpy saturated liquid.
 *
 * For pressure &gt; 16.5292 MPa the answer is obtained by iteration using
 * Ridders' root finding algorithm [Numerical Recipes, 3rd ed, 2007].
 *
 * @param pressure absolute pressure [MPa]
 * @return specific enthalpy [kJ/kg]
 */
func (r *Region) SpecificEnthalpySaturatedLiquidP(pressure float64) float64 {

	if pressure > ps13 {
		if pressure == constants.Pc {
			return constants.Hc
		}
		h := []float64{hs13, constants.Hc, math.NaN(), math.NaN()}
		dp := []float64{ps13 - pressure, constants.Pc - pressure, math.NaN(), math.NaN()}

		for i := 0; i < ITERATION_LIMIT; i++ {
			h[2] = (h[0] + h[1]) / 2.0
			dp[2] = SaturationPressureB34H(h[2]) - pressure
			h[3] = h[2] + (h[2]-h[0])*signum(dp[0]-dp[1])*dp[2]/math.Sqrt(dp[2]*dp[2]-dp[0]*dp[1])
			dp[3] = SaturationPressureB34H(h[3]) - pressure

			if dp[3] < 0 {
				h[0] = h[3]
				h[1] = h[2]
				dp[0] = dp[3]
				dp[1] = dp[2]
			} else {
				h[0] = h[2]
				h[1] = h[3]
				dp[0] = dp[2]
				dp[1] = dp[3]
			}
			if math.Abs(dp[3]) < TOLERANCE {
				break
			}
		}
		return h[3]

	} else {
		Ts := r.SaturationTemperatureP(pressure)

		return firstRegion.REGION1.SpecificEnthalpyPT(pressure, Ts)
	}
}

/**
 * Specific enthalpy saturated vapour.
 *
 * For pressure &gt; 16.5292 MPa the answer is obtained by iteration using
 * Ridders' root finding algorithm [Numerical Recipes, 3rd ed, 2007].
 *
 * @param pressure absolute pressure [MPa]
 * @return specific enthalpy [kJ/kg]
 */
func (r *Region) SpecificEnthalpySaturatedVapourP(pressure float64) float64 {

	if pressure > ps13 {
		if pressure == constants.Pc {
			return constants.Hc
		}
		h := []float64{constants.Hc, hs23, math.NaN(), math.NaN()}
		p := []float64{constants.Pc - pressure, ps13 - pressure, math.NaN(), math.NaN()}

		for i := 0; i < ITERATION_LIMIT; i++ {
			h[2] = (h[0] + h[1]) / 2.0
			p[2] = SaturationPressureB34H(h[2]) - pressure
			h[3] = h[2] + (h[2]-h[0])*signum(p[0]-p[1])*p[2]/math.Sqrt(p[2]*p[2]-p[0]*p[1])
			p[3] = SaturationPressureB34H(h[3]) - pressure

			if p[3] < 0 {
				h[0] = h[3]
				h[1] = h[2]
				p[0] = p[3]
				p[1] = p[2]
			} else {
				h[0] = h[2]
				h[1] = h[3]
				p[0] = p[2]
				p[1] = p[3]
			}
			if math.Abs(p[3]) < TOLERANCE {
				break
			}
		}
		return h[3]

	} else {
		Ts := r.SaturationTemperatureP(pressure)

		return secondRegion.REGION2.SpecificEnthalpyPT(pressure, Ts)
	}
}

func specificEntropiesP(pressure float64) []float64 {
	return []float64{
		REGION4.SpecificEntropySaturatedLiquidP(pressure),
		REGION4.SpecificEntropySaturatedVapourP(pressure),
	}
}

/**
 * Specific entropy as a function of pressure & specific enthalpy.
 *
 * @param pressure pressure [MPa]
 * @param enthalpy specific enthalpy [kJ/kg]
 * @return specific entropy [kJ/(kg K)]
 */
func (r *Region) SpecificEntropyPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.SpecificEntropyPX(pressure, x)
}

func (r *Region) SpecificEntropyPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

/**
 * Specific entropy as a function of pressure & vapour fraction.
 *
 * @param pressure pressure [MPa]
 * @param vapourFraction vapour fraction [-]
 * @return specific entropy [kJ/(kg K)]
 */
func (r *Region) SpecificEntropyPX(pressure float64, vapourFraction float64) float64 {

	s := specificEntropiesP(pressure)

	return valueX(vapourFraction, s)
}

func (r *Region) SpecificEntropyRhoT(density float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) SpecificEntropySaturatedLiquidP(pressure float64) float64 {

	Ts := r.SaturationTemperatureP(pressure)

	if pressure > ps13 {
		v := r.SpecificVolumeSaturatedLiquidP(pressure)

		return thirdRegion.REGION3.SpecificEntropyRhoT(1.0/v, Ts)
	}
	return firstRegion.REGION1.SpecificEntropyPT(pressure, Ts)
}

func (r *Region) SpecificEntropySaturatedVapourP(pressure float64) float64 {

	Ts := r.SaturationTemperatureP(pressure)

	if pressure > ps13 {
		v := r.SpecificVolumeSaturatedVapourP(pressure)

		return thirdRegion.REGION3.SpecificEntropyRhoT(1.0/v, Ts)
	}
	return secondRegion.REGION2.SpecificEntropyPT(pressure, Ts)
}

func (r *Region) SpecificGibbsFreeEnergyPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func specificInternalEnergiesP(pressure float64) []float64 {
	return []float64{
		REGION4.SpecificInternalEnergySaturatedLiquidP(pressure),
		REGION4.SpecificInternalEnergySaturatedVapourP(pressure),
	}
}

func (r *Region) SpecificInternalEnergyPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.SpecificInternalEnergyPX(pressure, x)
}

func (r *Region) SpecificInternalEnergyPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) SpecificInternalEnergyPX(pressure float64, vapourFraction float64) float64 {

	v := specificInternalEnergiesP(pressure)

	return valueX(vapourFraction, v)
}

func (r *Region) SpecificInternalEnergySaturatedLiquidP(pressure float64) float64 {

	Ts := r.SaturationTemperatureP(pressure)

	if pressure > ps13 {
		v := r.SpecificVolumeSaturatedLiquidP(pressure)

		return thirdRegion.REGION3.SpecificInternalEnergyRhoT(1.0/v, Ts)
	}
	return firstRegion.REGION1.SpecificInternalEnergyPT(pressure, Ts)
}

func (r *Region) SpecificInternalEnergySaturatedVapourP(pressure float64) float64 {

	Ts := r.SaturationTemperatureP(pressure)

	if pressure > ps13 {
		v := r.SpecificVolumeSaturatedVapourP(pressure)

		return thirdRegion.REGION3.SpecificInternalEnergyRhoT(1.0/v, Ts)
	}
	return secondRegion.REGION2.SpecificInternalEnergyPT(pressure, Ts)
}

func (r *Region) SpecificIsobaricHeatCapacityPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.SpecificIsobaricHeatCapacityPX(pressure, x)
}

func (r *Region) SpecificIsobaricHeatCapacityPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) SpecificIsobaricHeatCapacityPX(pressure float64, vapourFraction float64) float64 {

	T := r.SaturationTemperatureP(pressure)
	var cp []float64

	if pressure > ps13 {
		/*
		   Region 3
		*/
		h := specificEnthalpiesP(pressure)
		rho := densitiesRegion3(pressure, h)
		cp[0] = thirdRegion.REGION3.SpecificIsobaricHeatCapacityRhoT(rho[0], T)
		cp[1] = thirdRegion.REGION3.SpecificIsobaricHeatCapacityRhoT(rho[1], T)

	} else {
		/*
		   Regions 1 & 2
		*/
		cp[0] = firstRegion.REGION1.SpecificIsobaricHeatCapacityPT(pressure, T)
		cp[1] = secondRegion.REGION2.SpecificIsobaricHeatCapacityPT(pressure, T)
	}
	return valueX(vapourFraction, cp)
}

func (r *Region) SpecificIsochoricHeatCapacityPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.SpecificIsochoricHeatCapacityPX(pressure, x)
}

func (r *Region) SpecificIsochoricHeatCapacityPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) SpecificIsochoricHeatCapacityPX(pressure float64, vapourFraction float64) float64 {

	T := r.SaturationTemperatureP(pressure)
	var cv []float64

	if pressure > ps13 {
		/*
		   Region 3
		*/
		h := specificEnthalpiesP(pressure)
		rho := densitiesRegion3(pressure, h)
		cv[0] = thirdRegion.REGION3.SpecificIsochoricHeatCapacityRhoT(rho[0], T)
		cv[1] = thirdRegion.REGION3.SpecificIsochoricHeatCapacityRhoT(rho[1], T)

	} else {
		/*
		   Regions 1 & 2
		*/
		cv[0] = firstRegion.REGION1.SpecificIsochoricHeatCapacityPT(pressure, T)
		cv[1] = secondRegion.REGION2.SpecificIsochoricHeatCapacityPT(pressure, T)
	}
	return valueX(vapourFraction, cv)
}

func (r *Region) SpecificVolumeHS(enthalpy float64, entropy float64) float64 {

	p := r.PressureHS(enthalpy, entropy)
	x := r.VapourFractionHS(enthalpy, entropy)

	return r.SpecificVolumePX(p, x)
}

func (r *Region) SpecificVolumePH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.SpecificVolumePX(pressure, x)
}

func (r *Region) SpecificVolumePS(pressure float64, entropy float64) float64 {

	x := r.VapourFractionPS(pressure, entropy)

	return r.SpecificVolumePX(pressure, x)
}

func (r *Region) SpecificVolumePT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) SpecificVolumePX(pressure float64, vapourFraction float64) float64 {

	v := specificVolumesP(pressure)

	return valueX(vapourFraction, v)
}

/**
 * Specific volume saturated liquid.
 *
 * For pressure &gt; 16.5292 MPa the answer is obtained by iteration using
 * Van Wijngaarden/Dekker/Brent root finding algorithm [Numerical Recipes,
 * 3rd ed, 2007].
 *
 * @param pressure absolute pressure [MPa]
 * @return specific volume [m&sup3;/kg]
 */
func (r *Region) SpecificVolumeSaturatedLiquidP(pressure float64) float64 {

	Ts := r.SaturationTemperatureP(pressure)

	if pressure > ps13 {
		if pressure == constants.Pc {
			return 1.0 / constants.Rhoc
		}
		v := []float64{firstRegion.REGION1.SpecificVolumePT(ps13, constants.T13),
			1.0 / constants.Rhoc, math.NaN(), math.NaN(), math.NaN()}
		dp := []float64{thirdRegion.REGION3.PressureRhoT(1.0/v[0], Ts) - pressure,
			thirdRegion.REGION3.PressureRhoT(constants.Rhoc, Ts) - pressure, math.NaN()}

		v[2] = v[1]
		dp[2] = dp[1]

		for i := 0; i < ITERATION_LIMIT; i++ {
			if dp[1]*dp[2] > 0 {
				v[2] = v[0]
				dp[2] = dp[0]
				v[3] = v[1] - v[0]
				v[4] = v[1] - v[0]
			}
			if math.Abs(dp[2]) < math.Abs(dp[1]) {
				v[0] = v[1]
				v[1] = v[2]
				v[2] = v[0]
				dp[0] = dp[1]
				dp[1] = dp[2]
				dp[2] = dp[0]
			}
			tol1 := 2*math.SmallestNonzeroFloat64*math.Abs(v[1]) + TOLERANCE/2.0
			xm := (v[2] - v[1]) / 2.0

			if math.Abs(xm) <= tol1 || dp[1] == 0 {
				break
			}
			if math.Abs(v[4]) >= tol1 && math.Abs(dp[0]) > math.Abs(dp[1]) {
				s := dp[1] / dp[0]
				var p, q, r float64

				if v[0] == v[2] {
					p = 2 * xm * s
					q = 1 - s
				} else {
					q = dp[0] / dp[2]
					r = dp[1] / dp[2]
					p = s * (2*xm*q*(q-r) - (v[1]-v[0])*(r-1))
					q = (q - 1) * (r - 1) * (s - 1)
				}
				if p > 0 {
					q = -q
				}
				p = math.Abs(p)
				min1 := 3*xm*q - math.Abs(tol1*q)
				min2 := math.Abs(v[4] * q)

				if 2*p < math.Min(min1, min2) {
					v[4] = v[3]
					v[3] = p / q
				} else {
					v[3] = xm
					v[4] = v[3]
				}
			} else {
				v[3] = xm
				v[4] = v[3]
			}
			v[0] = v[1]
			dp[0] = dp[1]
			v[1] += func() float64 {
				if math.Abs(v[3]) > tol1 {
					return v[3]
				} else {
					return sign(tol1, xm)
				}

			}()

			dp[1] = thirdRegion.REGION3.PressureRhoT(1.0/v[1], Ts) - pressure
		}
		return v[1]

	} else {
		return firstRegion.REGION1.SpecificVolumePT(pressure, Ts)
	}
}

/**
 * Specific volume saturated vapour.
 *
 * For pressure &gt; 16.5292 MPa the answer is obtained by iteration using
 * Van Wijngaarden/Dekker/Brent root finding algorithm [Numerical Recipes,
 * 3rd ed, 2007].
 *
 * @param pressure absolute pressure [MPa]
 * @return specific volume [m&sup3;/kg]
 */
func (r *Region) SpecificVolumeSaturatedVapourP(pressure float64) float64 {

	Ts := r.SaturationTemperatureP(pressure)

	if pressure > ps13 {
		if pressure == constants.Pc {
			return 1.0 / constants.Rhoc
		}
		v := []float64{math.NaN(), secondRegion.REGION2.SpecificVolumePT(ps13, constants.T13), math.NaN(), math.NaN(), math.NaN()}
		dp := []float64{math.NaN(), thirdRegion.REGION3.PressureRhoT(1.0/v[1], constants.T13) - pressure, math.NaN()}

		/*
		   Bracket Root
		*/
		for i := 0; i < 1000; i++ {
			v[0] = v[1] - float64(i)*0.001
			dp[0] = thirdRegion.REGION3.PressureRhoT(1.0/v[0], Ts) - pressure

			if dp[0]*dp[1] < 0 {
				break
			}
		}
		v[2] = v[1]
		dp[2] = dp[1]

		for i := 0; i < ITERATION_LIMIT; i++ {
			if dp[1]*dp[2] > 0 {
				v[2] = v[0]
				dp[2] = dp[0]
				v[4], v[3] = v[1]-v[0], v[1]-v[0]
			}
			if math.Abs(dp[2]) < math.Abs(dp[1]) {
				v[0] = v[1]
				v[1] = v[2]
				v[2] = v[0]
				dp[0] = dp[1]
				dp[1] = dp[2]
				dp[2] = dp[0]
			}
			tol1 := 2*math.SmallestNonzeroFloat64*math.Abs(v[1]) + TOLERANCE/2.0
			xm := (v[2] - v[1]) / 2.0

			if math.Abs(xm) <= tol1 || dp[1] == 0 {
				break
			}
			if math.Abs(v[4]) >= tol1 && math.Abs(dp[0]) > math.Abs(dp[1]) {
				s := dp[1] / dp[0]
				var p, q, r float64

				if v[0] == v[2] {
					p = 2 * xm * s
					q = 1 - s
				} else {
					q = dp[0] / dp[2]
					r = dp[1] / dp[2]
					p = s * (2*xm*q*(q-r) - (v[1]-v[0])*(r-1))
					q = (q - 1) * (r - 1) * (s - 1)
				}
				if p > 0 {
					q = -q
				}
				p = math.Abs(p)
				min1 := 3*xm*q - math.Abs(tol1*q)
				min2 := math.Abs(v[4] * q)

				if 2*p < math.Min(min1, min2) {
					v[4] = v[3]
					v[3] = p / q
				} else {
					v[3] = xm
					v[4] = v[3]
				}
			} else {
				v[3] = xm
				v[4] = v[3]
			}
			v[0] = v[1]
			dp[0] = dp[1]
			v[1] += func() float64 {
				if math.Abs(v[3]) > tol1 {
					return v[3]
				} else {
					return sign(tol1, xm)
				}

			}()
			dp[1] = thirdRegion.REGION3.PressureRhoT(1.0/v[1], Ts) - pressure
		}
		return v[1]

	} else {
		return secondRegion.REGION2.SpecificVolumePT(pressure, Ts)
	}
}

func specificVolumesP(pressure float64) []float64 {
	return []float64{
		REGION4.SpecificVolumeSaturatedLiquidP(pressure),
		REGION4.SpecificVolumeSaturatedVapourP(pressure)}
}

func (r *Region) SpeedOfSoundPH(pressure float64, enthalpy float64) float64 {

	x := r.VapourFractionPH(pressure, enthalpy)

	return r.SpeedOfSoundPX(pressure, x)
}

func (r *Region) SpeedOfSoundPT(pressure float64, temperature float64) float64 {
	return math.NaN()
}

func (r *Region) SpeedOfSoundPX(pressure float64, vapourFraction float64) float64 {

	T := r.SaturationTemperatureP(pressure)
	h := specificEnthalpiesP(pressure)
	var w []float64

	if pressure > ps13 {
		/*
		   Region 3
		*/
		rho := densitiesRegion3(pressure, h)
		w[0] = thirdRegion.REGION3.SpeedOfSoundRhoT(rho[0], T)
		w[1] = thirdRegion.REGION3.SpeedOfSoundRhoT(rho[1], T)

	} else {
		/*
		   Regions 1 & 2
		*/
		w[0] = firstRegion.REGION1.SpeedOfSoundPT(pressure, T)
		w[1] = secondRegion.REGION2.SpeedOfSoundPT(pressure, T)
	}
	return valueX(vapourFraction, w)
}

/**
 * Surface tension as a function of temperature.
 *
 * @param temperature temperature [K]
 * @return surface tension [N/m]
 */
func (r *Region) SurfaceTensionT(temperature float64) float64 {

	theta := temperature / constants.Tc

	return 235.8 * math.Pow(1-theta, 1.256) * (1 - 0.625*(1-theta)) * 1e-3
}

/**
 * [IF97 Supplementary Release S04, June 2014]
 *
 * @param enthalpy specific enthalpy [kJ/kg]
 * @param entropy specific entropy [kJ/(kg K)]
 * @return saturation temperature [K]
 */
func (r *Region) TemperatureHS(enthalpy float64, entropy float64) float64 {

	var theta float64 = 0
	eta := enthalpy / 2800.0
	sigma := entropy / 9.2
	x := []float64{eta - 0.119, sigma - 1.07}

	for _, ijn := range IJnHS {
		theta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return theta * 550
}

func (r *Region) TemperaturePH(pressure float64, dummy float64) float64 {
	return r.SaturationTemperatureP(pressure)
}

func (r *Region) TemperaturePS(pressure float64, dummy float64) float64 {
	return r.SaturationTemperatureP(pressure)
}

func valueX(x float64, lim []float64) float64 {
	return lim[0] + x*(lim[1]-lim[0])
}

func vapourFraction(value float64, lim []float64) float64 {
	return (value - lim[0]) / (lim[1] - lim[0])
}

/**
 * [IF97 Supplementary Release S04, June 2014]
 *
 * Note that region 3 is not used.
 *
 * @param enthalpy specific enthalpy [kJ/kg]
 * @param entropy specific entropy [kJ/(kg K)]
 * @return vapour fraction [-]
 */
func (r *Region) VapourFractionHS(enthalpy float64, entropy float64) float64 {

	Ts := r.TemperatureHS(enthalpy, entropy)
	ps := r.SaturationPressureT(Ts)

	if ps > constants.Pc {
		return math.NaN()
	}
	h := []float64{
		secondRegion.REGION2.SpecificEnthalpyPT(ps, Ts),
		secondRegion.REGION2.SpecificEnthalpyPT(ps, Ts),
	}

	return vapourFraction(enthalpy, h)
}

func (r *Region) VapourFractionPH(pressure float64, enthalpy float64) float64 {

	if pressure > constants.Pc {
		return math.NaN()
	}
	h := specificEnthalpiesP(pressure)

	return vapourFraction(enthalpy, h)
}

func (r *Region) VapourFractionPS(pressure float64, entropy float64) float64 {

	if pressure > constants.Pc {
		return math.NaN()
	}
	s := specificEntropiesP(pressure)

	return vapourFraction(entropy, s)
}

func (r *Region) VapourFractionTS(temperature float64, entropy float64) float64 {

	ps := r.SaturationPressureT(temperature)

	if ps > constants.Pc {
		return math.NaN()
	}
	return r.VapourFractionPS(ps, entropy)
}
