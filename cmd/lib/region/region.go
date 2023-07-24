package region

import (
	"math"

	"if97.com/cmd/lib/fifthRegion"
	"if97.com/cmd/lib/firstRegion"
	"if97.com/cmd/lib/fourthRegion"
	rangeError "if97.com/cmd/lib/region/RangeError"
	"if97.com/cmd/lib/secondRegion"
	"if97.com/cmd/lib/thirdRegion"
	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
	"if97.com/cmd/lib/utils/units"
)

const (
	p5   = 50                  // upper pressure boundary of region 5 [MPa]
	p132 = 100                 // upper pressure boundary of regions 1, 3, and 2 [MPa]
	T13  = constants.T0 + 350  // temperature boundary between region 1 and 3 (623.15 K) [K]
	T25  = constants.T0 + 800  // temperature boundary between region 2 and 5 (1073.15 K) [K]
	T5   = constants.T0 + 2000 // upper temperature boundary of region 5 (2273.15 K) [K]
	
	s2bc = 5.85
)

var (
	
	p0   = fourthRegion.REGION4.SaturationPressureT(constants.T0)
	ps13 = fourthRegion.REGION4.SaturationPressureT(T13) // (16.529 MPa) [MPa]
	hs13 = firstRegion.REGION1.SpecificEnthalpyPT(ps13, T13)
	ss13 = firstRegion.REGION1.SpecificEntropyPT(ps13, T13)
	hs23 = secondRegion.REGION2.SpecificEnthalpyPT(ps13, T13)
	ss23 = secondRegion.REGION2.SpecificEntropyPT(ps13, T13)
	s2   = secondRegion.REGION2.SpecificEntropyPT(p0, constants.T0)
	sc = thirdRegion.REGION3.SpecificEntropyRhoT(constants.Rhoc, constants.Tc)
	nB23 = []float64{
		0.34805185628969e3,
		-.11671859879975e1,
		0.10192970039326e-2,
		0.57254459862746e3,
		0.13918839778870e2,
	}
)



type Region interface {
	GetName() string
	/**
	 * Heat capacity ratio.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return heat capacity ratio [-]
	 */
	HeatCapacityRatioPT(p float64, T float64) float64

	/**
	 * Isentropic exponent.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return isentropic exponent [-]
	 */
	IsentropicExponentPT(p float64, T float64) float64

	/**
	 * Isobaric cubic expansion coefficient.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return isobaric cubic expansion coefficient [1/K]
	 */
	IsobaricCubicExpansionCoefficientPT(p float64, T float64) float64
	/**
	 * Isothermal compressibility.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return isothermal compressibility [1/MPa]
	 */
	IsothermalCompressibilityPT(p float64, T float64) float64
	/**
	 * Pressure as a function of specific enthalpy & specific entropy.
	 *
	 * @param h specific enthalpy [kJ/kg]
	 * @param s specific entropy [kJ/kg-K]
	 * @return pressure [MPa]
	 */
	PressureHS(h float64, s float64) float64
	/**
	 * Specific enthalpy.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific enthalpy [kJ/kg]
	 */
	SpecificEnthalpyPT(p float64, T float64) float64

	/**
	 * Specific entropy.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific entropy [kJ/kg-K]
	 */
	SpecificEntropyPT(p float64, T float64) float64

	/**
	 * Specific entropy.
	 *
	 * @param rho density [kg/m&sup3;]
	 * @param T temperature [K]
	 * @return specific entropy [kJ/kg-K]
	 */
	SpecificEntropyRhoT(rho float64, T float64)float64

	/**
	 * Specific Gibbs Free Energy.
	 *
	 * @param  p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific entropy [kJ/kg]
	 */

	SpecificGibbsFreeEnergyPT(pressure float64, temperature float64) float64

	/**
	 * Specific internal energy.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific internal energy [kJ/kg]
	 */
	SpecificInternalEnergyPT(p float64, T float64) float64

	/**
	 * Specific isobaric heat capacity.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific isobaric heat capacity [kJ/kg-K]
	 */
	SpecificIsobaricHeatCapacityPT(p float64, T float64) float64

	/**
	 * Specific isochoric heat capacity.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific isochoric heat capacity [kJ/kg-K]
	 */
	SpecificIsochoricHeatCapacityPT(p float64, T float64) float64
	/**
	 * Specific volume.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return specific volume [m&sup3;/kg]
	 */
	SpecificVolumePT(p float64, T float64) float64
	/**
	* Specific volume.
	*
	* @param p pressure [MPa]
	* @param h specific enthalpy [kJ/kg]
	* @return specific volume [m&sup3;/kg]
	*/
	SpecificVolumePH(p float64, h float64) float64

	/**
	* Specific volume as a function of pressure & specific entropy.
	*
	* @param p pressure [MPa]
	* @param s specific entropy [kJ/kg-K]
	* @return specific volume [m&sup3;/kg]
	*/
	SpecificVolumePS(p float64, s float64) float64
	
	
	/**
	 * Speed of sound.
	 *
	 * @param p pressure [MPa]
	 * @param T temperature [K]
	 * @return speed of sound [m/s]
	 */
	SpeedOfSoundPT(p float64, T float64) float64
	/**
	 * Temperature.
	 *
	 * @param h specific enthalpy [kJ/kg]
	 * @param s specific entropy [kJ/kg-K]
	 * @return temperature [K]
	 */
	TemperatureHS(h float64, s float64) float64

	/**
	 * Temperature.
	 *
	 * @param p pressure [MPa]
	 * @param h specific enthalpy [kJ/kg]
	 * @return temperature [K]
	 */
	TemperaturePH(p float64, h float64) float64

	/**
	 * Temperature as a function of pressure & specific entropy.
	 *
	 * @param p pressconstants
	 * @param s specific entropy [kJ/kg-K]
	 * @return temperature [K]
	 */
	TemperaturePS(p float64, s float64) float64
}


var Select map[int]Region = map[int]Region{
	1:&firstRegion.REGION1,
	2:&secondRegion.REGION2,
	3:&thirdRegion.REGION3,
	4:&fourthRegion.REGION4,
	5:&fifthRegion.REGION5,	
}

type IF97Region struct {
	Name string;
	Region
}

var IF97reg IF97Region = IF97Region{
	"",
	&firstRegion.REGION1,
}








/**
 * Auxiliary equation for the boundary between regions 2 and 3.
 *
 * @param p pressure [MPa]
 * @return temperature [K]
 */
func temperatureB23P(p float64) float64 {
	return nB23[3] + math.Sqrt((p-nB23[4])/nB23[2])
}

/**
 * Auxiliary equation for the boundary between regions 2 and 3. [IF97
 * Supplementary Release S04, June 2014]
 *
 * @param h specific enthalpy [kJ/kg]
 * @param s specific entropy [kJ/kg-K]
 * @return temperature [K]
 */
func temperatureB23HS(h float64, s float64) float64 {

	var theta float64 = 0
	eta := h / 3e3
	sigma := s / 5.3
	x := []float64{eta - 0.727, sigma - 0.864, 900}
	IJn := [][]float64{
		{-12, 10, .629096260829810e-3},
		{-10, 8, -.823453502583165e-3},
		{-8, 3, .515446951519474e-7},
		{-4, 4, -.117565945784945e1},
		{-3, 3, .348519684726192e1},
		{-2, -6, -.507837382408313e-11},
		{-2, 2, -.284637670005479e1},
		{-2, 3, -.236092263939673e1},
		{-2, 4, .601492324973779e1},
		{0, 0, .148039650824546e1},
		{1, -3, .360075182221907e-3},
		{1, -2, -.126700045009952e-1},
		{1, 10, -.122184332521413e7},
		{3, -2, .149276502463272},
		{3, -1, .698733471798484},
		{5, -5, -.252207040114321e-1},
		{6, -6, .147151930985213e-1},
		{6, -3, -.108618917681849e1},
		{8, -8, -.936875039816322e-3},
		{8, -2, .819877897570217e2},
		{8, -1, -.182041861521835e3},
		{12, -12, .261907376402688e-5},
		{12, -1, -.291626417025961e5},
		{14, -12, .140660774926165e-4},
		{14, 1, .783237062349385e7}}

	for _, ijn := range IJn {
		theta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return theta * x[2]
}

/**
 * Specific enthalpy as a function of pressure & specific entropy.
 *
 * @param p pressure [MPa]
 * @param s specific entropy [kJ/kg-K]
 * @return specific enthalpy[kJ/kg]
 * @throws OutOfRangeException out-of-range exception
 */
func (r *IF97Region) SpecificEnthalpyPS(p float64, s float64) float64 {

	T := r.TemperaturePS(p, s)

	return r.SpecificEnthalpyPT(p, T)
}

/**
 * Specific volume as a function of specific enthalpy & specific entropy.
 *
 * @param h specific enthalpy [kJ/kg]
 * @param s specific entropy [kJ/kg-K]
 * @return specific volume [m&sup3;/kg]
 * @throws OutOfRangeException out-of-range exception
 */
func (r *IF97Region) SpecificVolumeHS(h float64, s float64) float64 {

	p := r.PressureHS(h, s)
	T := r.TemperaturePS(p, s)
	return r.SpecificVolumePT(p, T)
}




func saturationPressure3(s float64) float64 {

	var pi float64 = 0
	sigma := s / 5.2
	x := []float64{sigma - 1.03, sigma - 0.699}
	IJn := [][]float64{
		{0, 0, 0.639767553612785},
		{1, 1, -0.129727445396014e2},
		{1, 32, -0.224595125848403e16},
		{4, 7, 0.177466741801846e7},
		{12, 4, 0.717079349571538e10},
		{12, 14, -0.378829107169011e18},
		{16, 36, -0.955586736431328e35},
		{24, 10, 0.187269814676188e24},
		{28, 0, 0.119254746466473e12},
		{32, 18, 0.110649277244882e37}}

	for _, ijn := range IJn {
		pi += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return pi * 22
}

func specificEnthalpy1(s float64) float64 {

	var eta float64 = 0
	sigma := s / 3.8
	x := []float64{sigma - 1.09, sigma + 0.366e-4, 1700}
	IJn := [][]float64{
		{0, 14, .332171191705237},
		{0, 36, .611217706323496e-3},
		{1, 3, -.882092478906822e1},
		{1, 16, -.455628192543250},
		{2, 0, -.263483840850452e-4},
		{2, 5, -.223949661148062e2},
		{3, 4, -.428398660164013e1},
		{3, 36, -.616679338856916},
		{4, 4, -.146823031104040e2},
		{4, 16, .284523138727299e3},
		{4, 24, -.113398503195444e3},
		{5, 18, .115671380760859e4},
		{5, 24, .395551267359325e3},
		{7, 1, -.154891257229285e1},
		{8, 4, .194486637751291e2},
		{12, 2, -.357915139457043e1},
		{12, 4, -.335369414148819e1},
		{14, 1, -.664426796332460},
		{14, 22, .323321885383934e5},
		{16, 10, .331766744667084e4},
		{20, 12, -.223501257931087e5},
		{20, 28, .573953875852936e7},
		{22, 8, .173226193407919e3},
		{24, 3, -.363968822121321e-1},
		{28, 0, .834596332878346e-6},
		{32, 6, .503611916682674e1},
		{32, 8, .655444787064505e2}}

	for _, ijn := range IJn {
		eta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return eta * x[2]
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

/**
 * Specific enthalpy saturated liquid between region 3 and 4.
 *
 * @param s specific entropy [kJ/kg-K]
 * @return specific enthalpy [kJ/kg]
 */
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

func specificEnthalpyB13(s float64) float64 {

	var eta float64 = 0
	sigma := s / 3.8
	x := []float64{sigma - 0.884, sigma - 0.864, 1700}
	IJn := [][]float64{
		{0, 0, .913965547600543},
		{1, -2, -.430944856041991e-4},
		{1, 2, .603235694765419e2},
		{3, -12, .117518273082168e-17},
		{5, -4, .220000904781292},
		{6, -3, -.690815545851641e2}}

	for _, ijn := range IJn {
		eta += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return eta * x[2]
}

/**
 * Auxiliary equation for the boundary between regions 2 and 3.
 *
 * @param T temperature [K]
 * @return pressure [MPa]
 */
func pressureB23(T float64) float64 {
	return nB23[0] + nB23[1]*T + nB23[2]*T*T
}

/**
 * Get region as a function of specific enthalpy & specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return region
 * @throws OutOfRangeException out-of-range exception
 */
func (r *IF97Region) GetRegionHS(enthalpy float64, entropy float64) (int, error) {

	hB23limits := []float64{2.563592004e3, 2.812942061e3}
	sB23limits := []float64{5.048096828, 5.260578707}

	/*
	   Outer boundary Checks
	*/
	s1 := firstRegion.REGION1.SpecificEntropyPT(p132, constants.T0)
	s2 := secondRegion.REGION2.SpecificEntropyPT(p132, T25)

	h1 := []float64{
		firstRegion.REGION1.SpecificEnthalpyPT(p0, constants.T0),
		firstRegion.REGION1.SpecificEnthalpyPT(p132, constants.T0)}

	if enthalpy < h1[0] {
		return -1, rangeError.ErrorFromValue(quantity.H, enthalpy, h1[0])

	} else if entropy < 4.7516100567e-4 {
		p1 := firstRegion.REGION1.PressureHS(enthalpy, entropy)

		if firstRegion.REGION1.TemperaturePH(p1, enthalpy)+0.024 < constants.T0 {
			return -1, rangeError.ErrorFromValue(quantity.S, entropy,firstRegion.REGION1.SpecificEntropyPT(p1, constants.T0))

		}
	}
	if s1 <= entropy && entropy <= s2 {
		if entropy <= firstRegion.REGION1.SpecificEntropyPT(p132, T13) {
			h1Lim := firstRegion.REGION1.SpecificEnthalpyPT(
				p132, firstRegion.REGION1.TemperaturePS(p132, entropy),
			)

			if enthalpy > h1Lim {
				return -1, rangeError.ErrorFromValue(quantity.H, enthalpy, h1Lim)

			}
		} else if entropy <= secondRegion.REGION2.SpecificEntropyPT(p132, 863.15) {
			rho := 1.0 / thirdRegion.REGION3.SpecificVolumePS(p132, entropy)
			T := thirdRegion.REGION3.TemperaturePS(p132, entropy)
			hLim := thirdRegion.REGION3.SpecificEnthalpyRhoT(rho, T)

			if enthalpy > hLim {
				return -1, rangeError.ErrorFromValue(quantity.H, enthalpy, hLim)

			}
		}
	}
	/*
	   Select Region
	*/
	if entropy <= 3.778281340 {
		// region 1, 3, or 4
		if enthalpy <= specificEnthalpy1(entropy) {
			r.Region = &fourthRegion.REGION4
			return 4, nil

		} else if enthalpy > specificEnthalpyB13(entropy) {
			r.Region = &thirdRegion.REGION3
			return 3, nil

		} else {
			r.Region = &firstRegion.REGION1
			return 1, nil
		}
	} else if entropy <= sc {
		// region 3 or 4
		if enthalpy > specificEnthalpy3a(entropy) {
			r.Region = &thirdRegion.REGION3
			return 3, nil
		} else {
			r.Region = &fourthRegion.REGION4
			return 4, nil
		}
	} else if entropy < s2bc {
		if enthalpy <= specificEnthalpy2c3b(entropy) {
			r.Region = &fourthRegion.REGION4
			return 4, nil

		} else if enthalpy <= hB23limits[0] || entropy <= sB23limits[0] {
			r.Region = &thirdRegion.REGION3
			return 3, nil

		} else if enthalpy >= hB23limits[1] || entropy >= sB23limits[1] {
			r.Region = &secondRegion.REGION2
			return 2, nil

		} else if hB23limits[0] < enthalpy &&
			enthalpy < hB23limits[1] &&
			sB23limits[0] < entropy &&
			entropy < sB23limits[1] {
			if secondRegion.REGION2.PressureHS(enthalpy, entropy) >
				pressureB23(temperatureB23HS(enthalpy, entropy)) {
				r.Region = &thirdRegion.REGION3
				return 3, nil
			} else {
				r.Region = &secondRegion.REGION2
				return 2, nil
			}
		}
	} else if entropy <= 9.155759395 {
		if enthalpy <= specificEnthalpy2ab(entropy) {
			r.Region = &fourthRegion.REGION4
			return 4, nil
		}
	}
	r.Region = &secondRegion.REGION2
	return 2, nil
}

func (r *IF97Region) GetRegionPH(pressure float64, enthalpy float64) (int, error) {

	/*
	   Checks
	*/
	if pressure < p0 {
		return -1, rangeError.ErrorFromValue(quantity.P, pressure, p0);


	} else if pressure > p132 {
		return -1,  rangeError.ErrorFromValue(quantity.P, pressure, p132);

	}
	h25 := secondRegion.REGION2.SpecificEnthalpyPT(pressure, T25)

	/*
	   Select Region
	*/
	if enthalpy > h25 {
		if pressure > p5 {
			return -1, rangeError.RangeError{
				QUANTITIES: []quantity.Quantity{quantity.P, quantity.H}, 
				VALUES: []float64{pressure, enthalpy},
				LIMITS: []float64{p5, h25},
				UNIT_SYSTEM: units.DEFAULT,
			}

		}

	}
	if pressure <= ps13 {
		// region 1, 4, or 2
		Ts := fourthRegion.REGION4.SaturationTemperatureP(pressure)

		if enthalpy < firstRegion.REGION1.SpecificEnthalpyPT(pressure, Ts) {
			r.Region = &firstRegion.REGION1
			return 1, nil

		} else if enthalpy > secondRegion.REGION2.SpecificEnthalpyPT(pressure, Ts) {
			r.Region = &secondRegion.REGION2
			return 2, nil

		} else {
			r.Region = &fourthRegion.REGION4
			return 4, nil
		}
	} else if hs13 <= enthalpy && enthalpy <= hs23 {
		// region 3 or 4
		if pressure > fourthRegion.SaturationPressureB34H(enthalpy)*(1-4.3e-6) {
			r.Region = &thirdRegion.REGION3
			return 3, nil
		} else {
			r.Region = &fourthRegion.REGION4
			return 4, nil
		}

	} else if enthalpy <= firstRegion.REGION1.SpecificEnthalpyPT(pressure, T13) {
		r.Region = &firstRegion.REGION1
		return 1, nil

	} else if enthalpy >= secondRegion.REGION2.SpecificEnthalpyPT(pressure, temperatureB23P(pressure)) {
		r.Region = &secondRegion.REGION2
		return 2, nil

	} else {
		r.Region = &thirdRegion.REGION3
		return 3, nil
	}
}

func ConvertToDefault(quantity []float64, value float64) float64 {
	return value*quantity[0] + quantity[1]
}

/**
 * Returns the appropriate region.
 *
 * Note that this method never returns region 4, so prevent when possible.
 *
 * @param unitSystem unit system
 * @param pressure pressure [MPa]
 * @param temperature temperature [K]
 * @return region
 * @throws OutOfRangeException
 */
func (r *IF97Region) GetRegionPTUnits(unitSystem units.UnitSystem, pressure float64, temperature float64) (int, error) {

	p := ConvertToDefault(unitSystem.PRESSURE, pressure)
	T := ConvertToDefault(unitSystem.TEMPERATURE, temperature)

	retval, err := r.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	return retval, err
}

/**
 * Returns the appropriate region.
 *
 * Note that this method never returns region 4, so prevent when possible.
 *
 * @param pressure pressure [MPa]
 * @param temperature temperature [K]
 * @return region
 * @throws OutOfRangeException
 */
func (r *IF97Region) GetRegionPT(pressure float64, temperature float64) (int, error) {

	/*
	   Checks
	*/
	//return -1, ErrorFromValue(quantity.P, pressure, p0)
	if (pressure <= 0) {
	    return -1,rangeError.ErrorFromValue(quantity.P, pressure, 0)

	} else if (pressure > p132) {
	    return -1,rangeError.ErrorFromValue(quantity.P, pressure, p132)

	} else if (temperature < constants.T0) {
	    return -1,rangeError.ErrorFromValue(quantity.T, temperature, constants.T0)

	} else if (temperature > T25 && pressure > p5) {
	    return -1,rangeError.RangeError{
	    QUANTITIES: []quantity.Quantity{
	     quantity.P,
	    quantity.T,
	    }, 
		VALUES: []float64{pressure, temperature}, 
		LIMITS: []float64{p5, T25},
		 UNIT_SYSTEM: units.DEFAULT,
		}
	} else if (temperature > T5) {
	    return -1, rangeError.ErrorFromValue(quantity.T,temperature,T5)
	}

	/*
	   Select Region
	*/
	if temperature > T25 {
		r.Region = &fifthRegion.REGION5
		return 5, nil

	} else if temperature > T13 {
		if pressure > pressureB23(temperature) {
			r.Region = &thirdRegion.REGION3
			return 3, nil
		}
	} else if pressure > fourthRegion.REGION4.SaturationPressureT(temperature) {
		r.Region = &firstRegion.REGION1
		return 1, nil
	}
	r.Region = &secondRegion.REGION2
	return 2, nil
}

// /**
//   - Get region as a function of pressure & specific entropy.
//     *
//   - @param pressure pressure [MPa]
//   - @param entropy specific entropy [kJ/(kg K)]
//   - @return region
//   - @throws OutOfRangeException out-of-range exception
//     */
func (r *IF97Region) GetRegionPS(pressure float64, entropy float64) (int,error) {

	/*
	   Checks
	*/
	s1 := firstRegion.REGION1.SpecificEntropyPT(pressure, constants.T0)
	s2 := secondRegion.REGION2.SpecificEntropyPT(pressure, T25)
	//TODO : new Error
	if (pressure < p0) {
	    return -1, rangeError.ErrorFromValue(quantity.P, pressure, p0)

	} else if (pressure > p132) {
	    return-1, rangeError.ErrorFromValue(quantity.P, pressure, p132);

	} else if (entropy < s1) {
	    return-1, rangeError.ErrorFromValue(quantity.S, entropy, s1);

	} else if (entropy > s2) {
	    return-1, rangeError.ErrorFromValue(quantity.S, entropy, s2)
	}

	/*
	   Select Region
	*/
	if pressure < ps13 {
		Tsat := fourthRegion.REGION4.SaturationTemperatureP(pressure)

		if entropy < firstRegion.REGION1.SpecificEntropyPT(pressure, Tsat) {
			r.Region = &firstRegion.REGION1
			return 1,nil

		} else if entropy > secondRegion.REGION2.SpecificEntropyPT(pressure, Tsat) {
			r.Region = &secondRegion.REGION2
			return 2,nil

		} else {
			r.Region = &fourthRegion.REGION4
			return 4,nil
		}
	} else if firstRegion.REGION1.SpecificEntropyPT(ps13, T13) <= entropy &&
		entropy <= secondRegion.REGION2.SpecificEntropyPT(ps13, T13) &&
		pressure < saturationPressure3(entropy) {
			r.Region = &fourthRegion.REGION4
		return 4,nil

	} else if entropy <= firstRegion.REGION1.SpecificEntropyPT(pressure, T13) {
		r.Region = &firstRegion.REGION1
		return 1,nil

	} else if entropy < secondRegion.REGION2.SpecificEntropyPT(pressure, temperatureB23P(pressure)) {
		r.Region = &thirdRegion.REGION3
		return 3,nil

	} else {
		r.Region = &secondRegion.REGION2
		return 2,nil
	}
}
