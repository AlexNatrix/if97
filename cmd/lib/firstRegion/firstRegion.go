package firstRegion

import (
	"math"

	"if97.com/cmd/lib/utils/constants"
)

/**
 * Region 1.
 */

const (
	Name = "Region 1"
	Tref = 1386
	pRef = 16.53
)

var (
	IJnPT = [][]float64{
		{0, -2, +0.14632971213167},
		{0, -1, -0.84548187169114},
		{0, 00, -0.37563603672040e1},
		{0, 01, +0.33855169168385e1},
		{0, 02, -0.95791963387872},
		{0, 03, +0.15772038513228},
		{0, 04, -0.16616417199501e-1},
		{0, 05, +0.81214629983568e-3},
		{1, -9, +0.28319080123804e-3},
		{1, -7, -0.60706301565874e-3},
		{1, -1, -0.18990068218419e-1},
		{1, 00, -0.32529748770505e-1},
		{1, 01, -0.21841717175414e-1},
		{1, 03, -0.52838357969930e-4},
		{2, -3, -0.47184321073267e-3},
		{2, 00, -0.30001780793026e-3},
		{2, 01, +0.47661393906987e-4},
		{2, 03, -0.44141845330846e-5},
		{2, 17, -0.72694996297594e-15},
		{3, -4, -0.31679644845054e-4},
		{3, 00, -0.28270797985312e-5},
		{3, 06, -0.85205128120103e-9},
		{4, -5, -0.22425281908000e-5},
		{4, -2, -0.65171222895601e-6},
		{4, 10, -0.14341729937924e-12},
		{5, -8, -0.40516996860117e-6},
		{8, -11, -0.12734301741641e-8},
		{8, -06, -0.17424871230634e-9},
		{21, -29, -0.68762131295531e-18},
		{23, -31, +0.14478307828521e-19},
		{29, -38, +0.26335781662795e-22},
		{30, -39, -0.11947622640071e-22},
		{31, -40, +0.18228094581404e-23},
		{32, -41, -0.93537087292458e-25}}
	IJnHS = [][]float64{
		{0, 0, -.691997014660582},
		{0, 1, -.183612548787560e2},
		{0, 2, -.928332409297335e1},
		{0, 4, .659639569909906e2},
		{0, 5, -.162060388912024e2},
		{0, 6, .450620017338667e3},
		{0, 8, .854680678224170e3},
		{0, 14, .607523214001162e4},
		{1, 0, .326487682621856e2},
		{1, 1, -.269408844582931e2},
		{1, 4, -.319947848334300e3},
		{1, 6, -.928354307043320e3},
		{2, 0, .303634537455249e2},
		{2, 1, -.650540422444146e2},
		{2, 10, -.430991316516130e4},
		{3, 4, -.747512324096068e3},
		{4, 1, .730000345529245e3},
		{4, 4, .114284032569021e4},
		{5, 0, -.436407041874559e3}}
	IJnPH = [][]float64{
		{0, 00, -.23872489924521e3},
		{0, 01, 0.40421188637945e3},
		{0, 02, 0.11349746881718e3},
		{0, 06, -.58457616048039e1},
		{0, 22, -.15285482413140e-3},
		{0, 32, -.10866707695377e-5},
		{1, 00, -.13391744872602e2},
		{1, 01, 0.43211039183559e2},
		{1, 02, -.54010067170506e2},
		{1, 03, 0.30535892203916e2},
		{1, 04, -.65964749423638e1},
		{1, 10, 0.93965400878363e-2},
		{1, 32, 0.11573647505340e-6},
		{2, 10, -.25858641282073e-4},
		{2, 32, -.40644363084799e-8},
		{3, 10, 0.66456186191635e-7},
		{3, 32, 0.80670734103027e-10},
		{4, 32, -.93477771213947e-12},
		{5, 32, 0.58265442020601e-14},
		{6, 32, -.15020185953503e-16}}
	IJnPS = [][]float64{
		{0, 0, 0.17478268058307e3},
		{0, 1, 0.34806930892873e2},
		{0, 2, 0.65292584978455e1},
		{0, 3, 0.33039981775489},
		{0, 11, -0.19281382923196e-6},
		{0, 31, -0.24909197244573e-22},
		{1, 0, -0.26107636489332},
		{1, 1, 0.22592965981586},
		{1, 2, -0.64256463395226e-1},
		{1, 3, 0.78876289270526e-2},
		{1, 12, 0.35672110607366e-9},
		{1, 31, 0.17332496994895e-23},
		{2, 0, 0.56608900654837e-3},
		{2, 1, -0.32635483139717e-3},
		{2, 2, 0.44778286690632e-4},
		{2, 9, -0.51322156908507e-9},
		{2, 31, -0.42522657042207e-25},
		{3, 10, 0.26400441360689e-12},
		{3, 32, 0.78124600459723e-28},
		{4, 32, -0.30732199903668e-30}}
)

type Region struct {
	Name string
}

var REGION1 = Region{
	Name: Name,
}

func (r *Region)GetName() string{
	return r.Name
}

/**
 * Dimensionless Gibbs free energy.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return gamma
 */
func gamma(pi float64, tau float64) float64 {

	var out float64 = 0
	var x = []float64{7.1 - pi, tau - 1.222}

	for _, ijn := range IJnPT {
		out += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return out
}

/**
 * First partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return d/dpi gamma
 */
func gammaPi(pi float64, tau float64) float64 {

	var out float64 = 0
	var x = []float64{7.1 - pi, tau - 1.222}

	for _, ijn := range IJnPT {
		out -= ijn[2] * ijn[0] * math.Pow(x[0], ijn[0]-1) * math.Pow(x[1], ijn[1])
	}
	return out
}

/**
 * Second partial derivative with respect to pi.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return d2/dpi2 gamma
 */
func gammaPiPi(pi float64, tau float64) float64 {

	var out float64 = 0
	var x = []float64{7.1 - pi, tau - 1.222}

	for _, ijn := range IJnPT {
		out += ijn[2] * ijn[0] * (ijn[0] - 1) * math.Pow(x[0], ijn[0]-2) * math.Pow(x[1], ijn[1])
	}
	return out
}

/**
 * Second partial derivative with respect to pi & tau.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return d/dpi d/dtau gamma
 */
func gammaPiTau(pi float64, tau float64) float64 {

	var out float64 = 0
	var x = []float64{7.1 - pi, tau - 1.222}

	for _, ijn := range IJnPT {
		out -= ijn[2] * ijn[0] * math.Pow(x[0], ijn[0]-1) * ijn[1] * math.Pow(x[1], ijn[1]-1)
	}
	return out
}

/**
 * First partial derivative with respect to tau.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return d/dtau gamma
 */
func gammaTau(pi float64, tau float64) float64 {

	var out float64 = 0
	var x = []float64{7.1 - pi, tau - 1.222}

	for _, ijn := range IJnPT {
		out += ijn[2] * math.Pow(x[0], ijn[0]) * ijn[1] * math.Pow(x[1], ijn[1]-1)
	}
	return out
}

/**
 * Second partial derivative with respect to tau.
 *
 * @param pi dimensionless pressure
 * @param tau dimensionless temperature
 * @return d2/dtau2 gamma
 */
func gammaTauTau(pi float64, tau float64) float64 {

	var out float64 = 0
	var x = []float64{7.1 - pi, tau - 1.222}

	for _, ijn := range IJnPT {
		out += ijn[2] * math.Pow(x[0], ijn[0]) * ijn[1] * (ijn[1] - 1) * math.Pow(x[1], ijn[1]-2)
	}
	return out
}


func (r *Region)HeatCapacityRatioPT(pressure float64, temperature float64) float64 {

	var pi float64 = pressure / pRef
	var tau float64 = Tref / temperature
	var x float64 = gammaPi(pi, tau) - tau*gammaPiTau(pi, tau)

	return 1 / (1 - x*x/(tau*tau*gammaPiPi(pi, tau)*gammaTauTau(pi, tau)))
}

func (r *Region)IsentropicExponentPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	gPi := gammaPi(pi, tau)
	gPiPi := gammaPiPi(pi, tau)
	x := tau * tau * gammaTauTau(pi, tau)
	y := gPi - tau*gammaPiTau(pi, tau)

	return -gPi * x / (pi * (gPiPi*x - y*y))
}

func (r *Region)IsobaricCubicExpansionCoefficientPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (1 - tau*gammaPiTau(pi, tau)/gammaPi(pi, tau)) / temperature
}


func (r *Region)IsothermalCompressibilityPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return -pi * gammaPiPi(pi, tau) / gammaPi(pi, tau) / pressure
}

func (r *Region)PressureHS(enthalpy float64, entropy float64) float64 {

	var pi float64 = 0
	x := []float64{enthalpy/3400 + 0.05, entropy/7.6 + 0.05}

	for _, ijn := range IJnHS {
		pi += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1])
	}
	return pi * 100
}


func (r *Region)SpecificEnthalpyPT(pressure float64, temperature float64) float64 {

	tau := Tref / temperature

	return tau * gammaTau(pressure/pRef, tau) * constants.R * temperature
}


func (r *Region)SpecificEntropyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (tau*gammaTau(pi, tau) - gamma(pi, tau)) * constants.R
}


func (r *Region)SpecificEntropyRhoT(rho float64, T float64) float64 {
	panic("Region1.specificEntropyRhoT() pending implementation. Im sorry")
}


func (r *Region)SpecificGibbsFreeEnergyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return gamma(pi, tau) * constants.R * temperature
}


func (r *Region)SpecificInternalEnergyPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature

	return (tau*gammaTau(pi, tau) - pi*gammaPi(pi, tau)) * constants.R * temperature
}


func (r *Region)SpecificIsobaricHeatCapacityPT(pressure float64, temperature float64) float64 {

	tau := Tref / temperature

	return -tau * tau * gammaTauTau(pressure/pRef, tau) * constants.R
}


func (r *Region)SpecificIsochoricHeatCapacityPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	x := gammaPi(pi, tau) - tau*gammaPiTau(pi, tau)

	return (-tau*tau*gammaTauTau(pi, tau) + x*x/gammaPiPi(pi, tau)) * constants.R
}


func (r *Region)SpecificVolumePT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef

	return pi * gammaPi(pi, Tref/temperature) / 1e3 * constants.R * temperature / pressure
}


func (r *Region)SpeedOfSoundPT(pressure float64, temperature float64) float64 {

	pi := pressure / pRef
	tau := Tref / temperature
	gPi := gammaPi(pi, tau)
	x := gPi - tau*gammaPiTau(pi, tau)

	return math.Sqrt((gPi * gPi / (x*x/(tau*tau*gammaTauTau(pi, tau)) - gammaPiPi(pi, tau))) * 1e3 * constants.R * temperature)
}


func (r *Region)TemperatureHS(enthalpy float64, entropy float64) float64 {
	return r.TemperaturePH(r.PressureHS(enthalpy, entropy), enthalpy)
}


func (r *Region)TemperaturePH(pressure float64, enthalpy float64) float64 {

	var out float64 = 0
	x := enthalpy/2500 + 1

	for _, ijn := range IJnPH {
		out += ijn[2] * math.Pow(pressure, ijn[0]) * math.Pow(x, ijn[1])
	}
	return out
}


func (r *Region)TemperaturePS(pressure float64, entropy float64) float64 {

	var out float64 = 0

	for _, ijn := range IJnPS {
		out += ijn[2] * math.Pow(pressure, ijn[0]) * math.Pow(entropy+2, ijn[1])
	}
	return out
}


