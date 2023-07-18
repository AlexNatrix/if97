package calculator

import (
	"math"

	"if97.com/cmd/lib/fourthRegion"
	"if97.com/cmd/lib/thirdRegion"
	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
)

//FAKE!!!!!!!!!!!!!!
type Region interface {
	getRegionPH(p float64, h float64) sRegion
	specificVolumePH(p float64, h float64) float64
	temperaturePH(p float64, h float64) float64
	specificIsobaricHeatCapacityPT(p float64, T float64) float64
	specificVolumePT(p float64, T float64)float64
	specificEntropyPT(pMPa float64, T float64) float64
	isobaricCubicExpansionCoefficientPT(pMPa float64, T float64) float64
	isothermalCompressibilityPT(pMPa float64, T float64)float64

}

type sRegion struct {
	Name string
	Region
}

/**
 * Prandtl number.
 *
 * @param p pressure [MPa]
 * @param h specific enthalpy [kJ/kg]
 * @return Prandtl number [-]
 * @throws OutOfRangeException out-of-range exception
 */
func PrandtlPH(p float64, h float64) float64 {
	
	var region = Region.getRegionPH(p, h)

	var cp float64
	rho := 1 / region.specificVolumePH(p, h)
	T := region.temperaturePH(p, h)
	eta := DynamicViscosityRhoT(rho, T)
	lambda := ThermalConductivityRhoT(rho, T) / 1e3

	if region.Name == "Region4" {
		cp = fourthRegion.SpecificIsobaricHeatCapacityPH(p, h)

	} else {
		cp = region.specificIsobaricHeatCapacityPT(p, T)
	}
	return eta * cp / lambda
}

/**
 * Prandtl number.
 *
 * @param p pressure [MPa]
 * @param T temperature [K]
 * @return Prandtl number [-]
 * @throws OutOfRangeException out-of-range exception
 */
func PrandtlPT(p float64, T float64) float64 {

	var region = Region.getRegionPH(p, h)

	rho := 1 / region.specificVolumePT(p, T)
	eta := DynamicViscosityRhoT(rho, T)
	cp := region.specificIsobaricHeatCapacityPT(p, T)
	lambda := ThermalConductivityRhoT(rho, T) / 1e3

	return eta * cp / lambda
}

/**
 * Dielectric constant.
 *
 * @param rho density [kg/m3]
 * @param T temperature [K]
 * @return dielectric constant [-]
 * @throws OutOfRangeException out-of-range exception
 */
func DielectricConstantRhoT(rho float64, T float64) float64 {

	if T < 238.15 {
		panic("OOR") //OutOfRangeException(Quantity.T, T, 238.15);

	} else if T > 873.15 {
		panic("OOR") //OutOfRangeException(Quantity.T, T, 873.15);
	}

	k := 1.380658e-23
	NA := 6.0221367e23
	alpha := 1.636e-40
	epsilon0 := 8.854187817e-12
	mu := 6.138e-30
	M := 0.018015268
	n12 := 0.196096504426e-2
	delta := rho / constants.Rhoc
	tau := constants.Tc / T
	g := 1 + n12*delta*math.Pow(constants.Tc/228/tau-1, -1.2)
	IJn := [][]float64{
		{1, 0.25, 0.978224486826},
		{1, 1.0, -0.957771379375},
		{1, 2.5, 0.237511794148},
		{2, 1.5, 0.714692244396},
		{3, 1.5, -.298217036956},
		{3, 2.5, -.108863472196},
		{4, 2, 0.949327488264e-1},
		{5, 2, -.980469816509e-2},
		{6, 5, 0.165167634970e-4},
		{7, 0.5, 0.937359795772e-4},
		{10, 10, -.123179218720e-9}}

	for _, ijn := range IJn {
		g += ijn[2] * math.Pow(delta, ijn[0]) * math.Pow(tau, ijn[1])
	}

	A := NA * mu * mu * rho * g / (M * epsilon0 * k * T)
	B := NA * alpha * rho / (3 * M * epsilon0)

	return (1 + A + 5*B + math.Sqrt(9+2*A+18*B+A*A+10*A*B+9*B*B)) / (4 * (1 - B))
}

    /**
     * Dynamic viscosity.
     *
     * @param rho density [kg/m3]
     * @param T temperature [K]
     * @return dynamic viscosity [Pa-s]
     */
    func DynamicViscosityRhoT(rho float64, T float64) float64 {

        delta := rho / constants.Rhoc
                theta := T / constants.Tc
                 var psi0 float64 = 0
                var psi1 float64 = 0
         n0 := []float64{0.167752e-1, 0.220462e-1, 0.6366564e-2, -0.241605e-2};
         IJn := [][]float64{
            {0, 0, 0.520094},
            {0, 1, 0.850895e-1},
            {0, 2, -.108374e1},
            {0, 3, -.289555},
            {1, 0, 0.222531},
            {1, 1, 0.999115},
            {1, 2, 0.188797e1},
            {1, 3, 0.126613e1},
            {1, 5, 0.120573},
            {2, 0, -.281378},
            {2, 1, -.906851},
            {2, 2, -.772479},
            {2, 3, -.489837},
            {2, 4, -.257040},
            {3, 0, 0.161913},
            {3, 1, 0.257399},
            {4, 0, -.325372e-1},
            {4, 3, 0.698452e-1},
            {5, 4, 0.872102e-2},
            {6, 3, -.435673e-2},
            {6, 5, -.593264e-3}};

        for i := 0; i < len(n0); i++ {
            psi0 += n0[i] / math.Pow(theta, float64(i));
        }
        psi0 = math.Sqrt(theta) / psi0;

       x := []float64{delta - 1, 1 / theta - 1};

        for _,ijn := range IJn {
            psi1 += ijn[2] * math.Pow(x[0], ijn[0]) * math.Pow(x[1], ijn[1]);
        }
        psi1 = math.Exp(delta * psi1);

        return psi0 * psi1 * 1e-6;
    }

    /**
     * Gets the partial derivative of z with respect to x for constant y in
     * SI units, as a function of pressure and temperature in DEFAULT units.
     *
     * This method is for regions described by specific Gibbs free energy.
     *
     * @param region region
     * @param pMPa pressure [MPa]
     * @param T temperature [K]
     * @param x any quantity
     * @param y any quantity
     * @param z any quantity
     * @return partial derivative [SI units]
     * @throws OutOfRangeException out-of-range exception
     */
    func PartialDerivativePT(region Region,  pMPa float64 , T float64, x quantity.Quantity ,y quantity.Quantity, z quantity.Quantity) float64{

         		p := pMPa * 1e6 // [Pa]
                v := region.specificVolumePT(pMPa, T) // [m³/kg]
                s := region.specificEntropyPT(pMPa, T) * 1e3 // [J/(kg·K)]
                cp := region.specificIsobaricHeatCapacityPT(pMPa, T) * 1e3 // [J/(kg·K)]
                alphaV := region.isobaricCubicExpansionCoefficientPT(pMPa, T) // [1/K]
                kappaT := region.isothermalCompressibilityPT(pMPa, T) / 1e6; // [1/Pa]

        		dx := PartialDerivativesPT(p, T, x, v, s, cp, alphaV, kappaT)// [SI units]
                dy := PartialDerivativesPT(p, T, y, v, s, cp, alphaV, kappaT) // [SI units]
                dz := PartialDerivativesPT(p, T, z, v, s, cp, alphaV, kappaT); // [SI units]

        		 dx_dT := dx[0]
                dy_dT := dy[0]
                dz_dT := dz[0]
                dx_dp := dx[1]
                dy_dp := dy[1]
                dz_dp := dz[1];

        return (dz_dp * dy_dT - dz_dT * dy_dp) / (dx_dp * dy_dT - dx_dT * dy_dp);
    }

    /**
     * Gets the partial derivative of z with respect to x for constant y in
     * SI units, as a function of specific volume and temperature.
     *
     * This method is for region 3 described by specific Helmholtz free
     * energy.
     *
     * @param rho density [kg/m³]
     * @param T temperature [K]
     * @param x any quantity
     * @param y any quantity
     * @param z any quantity
     * @return partial derivative [SI units]
     */
    func PartialDerivativeRhoT( rho float64, T float64, x quantity.Quantity, y quantity.Quantity,  z quantity.Quantity) float64{

        		 v := 1 / rho // [m³/kg]
                p := thirdRegion.PressureRhoT(rho, T) * 1e6 // [Pa]
                s := thirdRegion.SpecificEntropyRhoT(rho, T) * 1e3 // [J/(kg·K)]]
                cv := thirdRegion.SpecificIsochoricHeatCapacityRhoT(rho, T) * 1e3 // [J/(kg·K)]
                alphap := thirdRegion.RelativePressureCoefficientRhoT(rho, T) // [1/K]
                betap := thirdRegion.IsothermalStressCoefficientRhoT(rho, T); // [kg/m³]

        		dx := PartialDerivativesVT(v, T, x, p, s, cv, alphap, betap)
                dy := PartialDerivativesVT(v, T, y, p, s, cv, alphap, betap)
                dz := PartialDerivativesVT(v, T, z, p, s, cv, alphap, betap)

        		 dx_dv := dx[0]
                dy_dv := dy[0]
                dz_dv := dz[0]
                dx_dT := dx[1]
                dy_dT := dy[1]
                dz_dT := dz[1];

        return (dz_dv * dy_dT - dz_dT * dy_dv) / (dx_dv * dy_dT - dx_dT * dy_dv);
    }

    /**
     * Partial derivatives of the given quantity with respect to pressure
     * and temperature for regions described by specific Gibbs free energy
     * in SI units.
     *
     * @param p pressure [Pa]
     * @param T temperature [K]
     * @param quantity quantity
     * @param v specific volume [m³/kg]
     * @param s specific entropy [J/(kg·K)]
     * @param cp specific isobaric heat capacity [J/(kg·K)]
     * @param alphaV isobaric cubic expansion coefficient [1/K]
     * @param kappaT isothermal compressibility [1/Pa]
     * @return partial derivatives d/dT and d/dp [SI units]
     */
    func PartialDerivativesPT( p float64,  T float64,  quantity quantity.Quantity,  v float64,  s float64,  cp float64, alphaV float64,  kappaT float64) []float64 {

        var d_dT, d_dp float64;

        switch (quantity.SYMBOL) {
            case "p":
                d_dT = 0; // [-]
                d_dp = 1; // [-]


			case "T":
                d_dT = 1; // [-]
                d_dp = 0; // [-]
 
            case "v":
                d_dT = v * alphaV; // [m³/(kg·K)]
                d_dp = -v * kappaT; // [m³/(kg·Pa)]

            case "u":
                d_dT = cp - p * v * alphaV; // [J/(kg·K)]
                d_dp = v * (p * kappaT - T * alphaV); // [m³/kg]
                break;

            case "h":
                d_dT = cp; // [J/(kg·K)]
                d_dp = v * (1 - T * alphaV); // [m³/kg]


            case "s":
                d_dT = cp / T; // [J/(kg·K²)]
                d_dp = -v * alphaV; // [m³/(kg·K)]


            case "g":
                d_dT = -s; // [J/(kg·K)]
                d_dp = v; // [m³/kg]
      

            case "f":
                d_dT = -p * v * alphaV - s; // [J/(kg·K)]
                d_dp = p * v * kappaT; // [m³/kg]
       

            case "\u03c1":
                d_dT = -alphaV / v; // [kg/(m³K)]
                d_dp = kappaT / v; // [kg/(m³Pa)]
             

            default:
                panic("Unsupported quantity for partial derivative: Calculator PartialDerivativePT" );
        }

        return []float64{d_dT, d_dp};
    }

    /**
     * Gets the partial derivatives of the given quantity with respect to
     * specific volume and temperature in SI units, for region 3 described
     * by specific Helmholtz free energy.
     *
     * @param v specific volume [m³/kg]
     * @param T temperature [K]
     * @param quantity quantity
     * @param p pressure [Pa]
     * @param s specific entropy [J/(kg·K)]
     * @param cv specific isochoric heat capacity [J/(kg·K)]
     * @param alphap relative pressure coefficient [1/K]
     * @param betap isothermal stress coefficient [kg/m³]
     * @return partial derivatives d/d&nu; and d/dT [SI units]
     */
    func PartialDerivativesVT( v float64, T float64, quantity quantity.Quantity, p float64, s float64, cv float64, alphap float64,  betap float64) []float64 {

        var d_dv, d_dT float64

        switch (quantity.SYMBOL) {
            case "p":
                d_dv = -p * betap; // [Pa·kg/m³]
                d_dT = p * alphap; // [Pa/K]
                break;

            case "T":
                d_dv = 0; // [-]
                d_dT = 1; // [-]
                break;

            case "v":
                d_dv = 1; // [-]
                d_dT = 0; // [-]
                break;

            case "u":
                d_dv = p * (T * alphap - 1); // [Pa]
                d_dT = cv; // [J/(kg·K)]
                break;

            case "h":
                d_dv = p * (T * alphap - v * betap); // [Pa]
                d_dT = cv + p * v * alphap; // [J/(kg·K)]
                break;

            case "s":
                d_dv = p * alphap; // [Pa/K]
                d_dT = cv / T; // [J/(kg·K²)]
                break;

            case "g":
                d_dv = -p * v * betap; // [Pa]
                d_dT = p * v * alphap - s; // [J/(kg·K)]
                break;

            case "f":
                d_dv = -p; // [Pa]
                d_dT = -s; // [J/(kg·K)]
                break;

            case "\u03c1":
                d_dv = -1 / (v * v); // [kg²/m^6]
                d_dT = 0; // [-]
                break;

            default:
                panic("Unsupported quantity for partial derivative: Calculator partialDerivativesVT" );
        }
        return []float64{d_dv, d_dT};
    }

    /**
     * Refractive index.
     *
     * @param rho density [kg/m&sup3;]
     * @param T temperature [K]
     * @param lambdaL wavelength [&mu;m]
     * @return refractive index [-]
     * @throws OutOfRangeException out-of-range exception
     */
    func  RefractiveIndexRhoTLambda(rho float64, T float64, lambdaL float64) float64 {

        if (T < 261.15) {
            panic("RefractiveIndexRhoTLambda") //OutOfRangeException(Quantity.T, T, 261.15);

        } else if (T > 773.15) {
            panic("RefractiveIndexRhoTLambda")  //OutOfRangeException(Quantity.T, T, 773.15);

        } else if (rho <= 0) {
            panic("RefractiveIndexRhoTLambda")  //OutOfRangeException(Quantity.rho, rho, 0);

        } else if (rho > 1060) {
            panic("RefractiveIndexRhoTLambda")  //OutOfRangeException(Quantity.rho, rho, 1060);

        } else if (lambdaL < 0.2) {
            panic("RefractiveIndexRhoTLambda")  //OutOfRangeException(Quantity.lambdaL, lambdaL, 0.2);

        } else if (lambdaL > 1.1) {
            panic("RefractiveIndexRhoTLambda")  //OutOfRangeException(Quantity.lambdaL, lambdaL, 1.1);
        }

         a := []float64{0.244257733, 0.974634476e-2, -.373234996e-2, 0.268678472e-3, 0.158920570e-2, 0.245934259e-2, 0.900704920, -.166626219e-1};

        	 delta := rho / 1e3
                theta := T / constants.T0
                Lambda := lambdaL / 0.589
                Lambda2 := Lambda * Lambda
                LambdaIR := 5.432937
                LambdaUV := -0.229202
                A := delta * (a[0] + a[1] * delta + a[2] * theta + a[3] * Lambda2 * theta + a[4] / Lambda2 + a[5] / (Lambda2 - LambdaUV * LambdaUV) + a[6] / (Lambda2 - LambdaIR * LambdaIR) + a[7] * delta * delta);

        return math.Sqrt((2 * A + 1) / (1 - A));
    }

    func ThermalConductivityRhoT(rho float64, T float64) float64 {

        /*
         Coefficients
         */
        	n0 := []float64{0.102811e-1, 0.299621e-1, 0.156146e-1, -.422464e-2}
                n1 := []float64{-.397070, 0.400302, 0.106000e1, -.171587, 0.239219e1}
                n2 := []float64{0.701309e-1, 0.118520e-1, 0.642857, 0.169937e-2, -.102000e1, -.411717e1, -.617937e1, 0.822994e-1, 0.100932e2, 0.308976e-2};

        		theta := T / 647.26
                DeltaTheta := math.Abs(theta - 1) + n2[9]
                delta := rho / 317.7
                var Lambda0 float64 = 0
				var A float64
                B := 2 + n2[7] * math.Pow(DeltaTheta, -0.6);

        for  i := 0; i < 4; i++ {
            Lambda0 += n0[i] * math.Pow(theta, float64(i));
        }
		if theta < 1 {
			A = n2[8]/ math.Pow(DeltaTheta, 0.6) 
		} else {
			A = 1 / DeltaTheta;
		}
          

        		Lambda1 := n1[0] + n1[1] * delta + n1[2] * math.Exp(n1[3] * math.Pow(delta + n1[4], 2))
                Lambda2 := (n2[0] / math.Pow(theta, 10) + n2[1]) * math.Pow(delta, 1.8) * math.Exp(n2[2] * (1.0 - math.Pow(delta, 2.8))) + 
				n2[3] * A * math.Pow(delta, B) * math.Exp(B / (1.0 + B) * (1.0 - math.Pow(delta, 1.0 + B))) + n2[4] * math.Exp(n2[5] * math.Pow(theta, 1.5) + n2[6] / math.Pow(delta, 5));

        return math.Sqrt(theta) * Lambda0 + Lambda1 + Lambda2;
    }

    func  ThermalDiffusivityPH(p float64, h float64) float64 {

         	var region  sRegion = Region.getRegionPH(p, h);

         		 rho := 1 / region.specificVolumePH(p, h)
                T := region.temperaturePH(p, h)
                lambda := ThermalConductivityRhoT(rho, T)
                var cp float64;

        if (region.Name == "Region4") {
            cp = fourthRegion.SpecificIsobaricHeatCapacityPH(p, h);

        } else {
            cp = region.specificIsobaricHeatCapacityPT(p, T);
        }
        return lambda / rho / cp;
    }
