package IF97

import (
	"errors"
	"math"

	"if97.com/cmd/lib/firstRegion"
	"if97.com/cmd/lib/fourthRegion"
	"if97.com/cmd/lib/region"
	rangeError "if97.com/cmd/lib/region/RangeError"
	"if97.com/cmd/lib/thirdRegion"
	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
	"if97.com/cmd/lib/utils/units"
)

type UNITSYSNUM int

const (
	_ UNITSYSNUM  = iota
	DEFAULT
)

type WSPcalculator interface {
}



type IF97 struct {
	UnitSystem units.UnitSystem;
	*region.IF97Region
}

var IF97test IF97calc = region.IF97Region{}

func New(unitSystem UNITSYSNUM) IF97{
	return IF97{
		units.UNITSYSTEMS[int(unitSystem)],
		&region.IF97Region{"",&firstRegion.REGION1},
	}
}

func (if97 *IF97) SetUnitSystem(unitsysnum UNITSYSNUM){
	if97.UnitSystem = units.UNITSYSTEMS[int(unitsysnum) ]
}

func (if97 *IF97)ConvertToDefault(quantity []float64, value float64) float64 {
	return value*quantity[0] + quantity[1]
}


func (if97 *IF97)ConvertFromDefault(quantity []float64, value float64)float64 {
	return (value - quantity[1]) / quantity[0];
}

func (if97 *IF97)ConvertFromDefaultQuantity(unitSystem units.UnitSystem, quantity quantity.Quantity, value float64) (float64, error) {

	switch quantity.NAME {
	case "temperature":
		return if97.ConvertFromDefault(unitSystem.TEMPERATURE, value),nil

	case "specific Helmholtz free energy":
	case "specific Gibbs free energy":
	case "specific internal energy":
		return if97.ConvertFromDefault(unitSystem.SPECIFIC_ENERGY, value),nil

	case "specific enthalpy":
		return if97.ConvertFromDefault(unitSystem.SPECIFIC_ENTHALPY, value),nil

	case "thermal conductivity":
		return if97.ConvertFromDefault(unitSystem.THERMAL_CONDUCTIVITY, value),nil

	case "wavelength":
		return if97.ConvertFromDefault(unitSystem.WAVELENGTH, value),nil

	case "absolute pressure":
		return if97.ConvertFromDefault(unitSystem.PRESSURE, value),nil

	case "density":
		return if97.ConvertFromDefault(unitSystem.DENSITY, value),nil

	case "specific entropy":
		return if97.ConvertFromDefault(unitSystem.SPECIFIC_ENTROPY, value),nil

	case "specific volume":
		return if97.ConvertFromDefault(unitSystem.SPECIFIC_VOLUME, value),nil

	case "vapour fraction":
		return value,nil

	}
	return -1, errors.New("no conversion available for: " + quantity.String())
}

/**
 * Prandtl number.
 *
 * @param p pressure [MPa]
 * @param specific entropy [kJ/K-kg]
 * @return Prandtl number [-]
 * @return  RangeError 
 */

func (if97 *IF97)PrandtlHS(enthalpy float64, entropy float64) (float64 , error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	_,err := if97.GetRegionHS(h, s);
	if err != nil{
		return -1,err
	}
	p := if97.Region.PressureHS(h, s)
	T := if97.Region.TemperatureHS(h, s)

	return if97.calcPrandtlPT(p, T);
}


/**
 * Prandtl number.
 *
 * @param p pressure [MPa]
 * @param h specific enthalpy [kJ/kg]
 * @return Prandtl number [-]
 * @return  RangeError 
 */
func (if97 *IF97)PrandtlPH(pressure float64, enthalpy float64) (float64, error) {

	 p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	 h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)

	return if97.calcPrandtlPH(p, h);
}




//helper function
 func (if97 *IF97) calcPrandtlPH(p float64, h float64) (float64,error) {
	regNum,err := if97.GetRegionPH(p, h)
    if err != nil {
        return -1,err
    }
	reg := region.IF97Region{
        Name:region.Select[regNum].GetName(),
        Region:region.Select[regNum],
    }

	var cp float64
	
	rho := 1 / reg.SpecificVolumePH(p, h)
	T := reg.TemperaturePH(p, h)
	eta,_ := if97.DynamicViscosityRhoT(rho, T)
	lambda := if97.ThermalConductivityRhoT(rho, T) / 1e3

	if reg.Name == "Region4" {
		cp = fourthRegion.REGION4.SpecificIsobaricHeatCapacityPH(p, h)

	} else {
		cp = reg.SpecificIsobaricHeatCapacityPT(p, T)
	}
	return eta * cp / lambda,nil
}

//helper function
func (if97 *IF97)calcPrandtlPT(p float64, T float64) (float64,error) {
	regNum,err := if97.GetRegionPT(p, T)
    if err != nil {
        return -1,err
    }
	reg := region.IF97Region{
        Name:region.Select[regNum].GetName(),
        Region:region.Select[regNum],
    }

	rho := 1 / reg.SpecificVolumePT(p, T)
	eta,_ := if97.DynamicViscosityRhoT(rho, T)
	cp := reg.SpecificIsobaricHeatCapacityPT(p, T)
	lambda := if97.ThermalConductivityRhoT(rho, T) / 1e3

	return eta * cp / lambda,nil
}

/**
 * Dielectric constant.
 *
 * @param rho density [kg/m3]
 * @param T temperature [K]
 * @return dielectric constant [-]
 * @return  RangeError 
 */
func (if97 *IF97)DielectricConstantRhoT(rho float64, T float64) (float64,error) {

	if T < 238.15 {
		return -1,rangeError.ErrorFromValue(quantity.T, T, 238.15);

	} else if T > 873.15 {
		return -1,rangeError.ErrorFromValue(quantity.T, T, 873.15);
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

	return (1 + A + 5*B + math.Sqrt(9+2*A+18*B+A*A+10*A*B+9*B*B)) / (4 * (1 - B)),nil
}

    /**
     * Dynamic viscosity.
     *
     * @param rho density [kg/m3]
     * @param T temperature [K]
     * @return dynamic viscosity [Pa-s]
     */
    func (if97 *IF97) DynamicViscosityRhoT(rho float64, T float64) (float64,error) {

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

        return psi0 * psi1 * 1e-6,nil
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
     * @return  RangeError 
     */
    func PartialDerivativePT(region region.IF97Region,  pMPa float64 , T float64, x quantity.Quantity ,y quantity.Quantity, z quantity.Quantity) (float64,error){

         		p := pMPa * 1e6 // [Pa]
                v := region.SpecificVolumePT(pMPa, T) // [m³/kg]
                s := region.SpecificEntropyPT(pMPa, T) * 1e3 // [J/(kg·K)]
                cp := region.SpecificIsobaricHeatCapacityPT(pMPa, T) * 1e3 // [J/(kg·K)]
                alphaV := region.IsobaricCubicExpansionCoefficientPT(pMPa, T) // [1/K]
                kappaT := region.IsothermalCompressibilityPT(pMPa, T) / 1e6; // [1/Pa]

        		dx,err := PartialDerivativesPT(p, T, x, v, s, cp, alphaV, kappaT)// [SI units]
				if err!=nil{
					return-1, err
				}
                dy,err := PartialDerivativesPT(p, T, y, v, s, cp, alphaV, kappaT) // [SI units]
				if err!=nil{
					return-1, err
				}
                dz,err := PartialDerivativesPT(p, T, z, v, s, cp, alphaV, kappaT); // [SI units]
				if err!=nil{
					return-1, err
				}

        		 dx_dT := dx[0]
                dy_dT := dy[0]
                dz_dT := dz[0]
                dx_dp := dx[1]
                dy_dp := dy[1]
                dz_dp := dz[1];

        return (dz_dp * dy_dT - dz_dT * dy_dp) / (dx_dp * dy_dT - dx_dT * dy_dp),nil;
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
    func PartialDerivativeRhoT( rho float64, T float64, x quantity.Quantity, y quantity.Quantity,  z quantity.Quantity) (float64,error){

        		 v := 1 / rho // [m³/kg]
                p := thirdRegion.REGION3.PressureRhoT(rho, T) * 1e6 // [Pa]
                s := thirdRegion.REGION3.SpecificEntropyRhoT(rho, T) * 1e3 // [J/(kg·K)]]
                cv := thirdRegion.REGION3.SpecificIsochoricHeatCapacityRhoT(rho, T) * 1e3 // [J/(kg·K)]
                alphap := thirdRegion.REGION3.RelativePressureCoefficientRhoT(rho, T) // [1/K]
                betap := thirdRegion.REGION3.IsothermalStressCoefficientRhoT(rho, T); // [kg/m³]

        		dx,err := PartialDerivativesVT(v, T, x, p, s, cv, alphap, betap)
				if err!=nil{
					return-1, err
				}
                dy,err := PartialDerivativesVT(v, T, y, p, s, cv, alphap, betap)
				if err!=nil{
					return-1, err
				}
                dz,err := PartialDerivativesVT(v, T, z, p, s, cv, alphap, betap)
				if err!=nil{
					return-1, err
				}

        		 dx_dv := dx[0]
                dy_dv := dy[0]
                dz_dv := dz[0]
                dx_dT := dx[1]
                dy_dT := dy[1]
                dz_dT := dz[1];

        return (dz_dv * dy_dT - dz_dT * dy_dv) / (dx_dv * dy_dT - dx_dT * dy_dv),nil;
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
    func PartialDerivativesPT( p float64,  T float64,  quantity quantity.Quantity,  v float64,  s float64,  cp float64, alphaV float64,  kappaT float64) ([]float64,error) {

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
                return []float64{-1, -1},errors.New("unsupported quantity for partial derivative: Calculator PartialDerivativePT" );
        }

        return []float64{d_dT, d_dp},nil;
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
    func PartialDerivativesVT( v float64, T float64, quantity quantity.Quantity, p float64, s float64, cv float64, alphap float64,  betap float64) ([]float64,error) {

        var d_dv, d_dT float64

        switch (quantity.SYMBOL) {
            case "p":
                d_dv = -p * betap; // [Pa·kg/m³]
                d_dT = p * alphap; // [Pa/K]
   

            case "T":
                d_dv = 0; // [-]
                d_dT = 1; // [-]
  

            case "v":
                d_dv = 1; // [-]
                d_dT = 0; // [-]


            case "u":
                d_dv = p * (T * alphap - 1); // [Pa]
                d_dT = cv; // [J/(kg·K)]


            case "h":
                d_dv = p * (T * alphap - v * betap); // [Pa]
                d_dT = cv + p * v * alphap; // [J/(kg·K)]
    

            case "s":
                d_dv = p * alphap; // [Pa/K]
                d_dT = cv / T; // [J/(kg·K²)]
    

            case "g":
                d_dv = -p * v * betap; // [Pa]
                d_dT = p * v * alphap - s; // [J/(kg·K)]

            case "f":
                d_dv = -p; // [Pa]
                d_dT = -s; // [J/(kg·K)]

            case "\u03c1":
                d_dv = -1 / (v * v); // [kg²/m^6]
                d_dT = 0; // [-]

            default:
                return []float64{-1, -1},errors.New("unsupported quantity for partial derivative: Calculator partialDerivativesVT" );
        }
        return []float64{d_dv, d_dT},nil;
    }

    /**
     * Refractive index.
     *
     * @param rho density [kg/m&sup3;]
     * @param T temperature [K]
     * @param lambdaL wavelength [&mu;m]
     * @return refractive index [-]
     * @return  RangeError 
     */
    func  (if97 *IF97)RefractiveIndexRhoTLambda(rho float64, T float64, lambdaL float64) (float64,error) {

        if (T < 261.15) {
            return -1,rangeError.ErrorFromValue(quantity.T, T, 261.15);

        } else if (T > 773.15) {
			return -1,rangeError.ErrorFromValue(quantity.T, T,  773.15)

        } else if (rho <= 0) {
            return -1,rangeError.ErrorFromValue(quantity.Rho, rho, 0) 

        } else if (rho > 1060) {
			return -1,rangeError.ErrorFromValue(quantity.Rho, rho, 1060)

        } else if (lambdaL < 0.2) {
			return -1,rangeError.ErrorFromValue(quantity.LambdaL, rho, 0.2)

        } else if (lambdaL > 1.1) {
			return -1,rangeError.ErrorFromValue(quantity.LambdaL, rho,  1.1)
        }

         a := []float64{0.244257733, 0.974634476e-2, -.373234996e-2, 0.268678472e-3, 0.158920570e-2, 0.245934259e-2, 0.900704920, -.166626219e-1};

        	 delta := rho / 1e3
                theta := T / constants.T0
                Lambda := lambdaL / 0.589
                Lambda2 := Lambda * Lambda
                LambdaIR := 5.432937
                LambdaUV := -0.229202
                A := delta * (a[0] + a[1] * delta + a[2] * theta + a[3] * Lambda2 * theta + a[4] / Lambda2 + a[5] / (Lambda2 - LambdaUV * LambdaUV) + a[6] / (Lambda2 - LambdaIR * LambdaIR) + a[7] * delta * delta);

        return math.Sqrt((2 * A + 1) / (1 - A)),nil;
    }


    /**
     * Dielectric constant as a function of pressure and temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return dielectric constant [-]
     * @throws OutOfRangeException out-of-range exception
     * @see #dielectricConstantRhoT(double, double)
     */
	 func (if97 *IF97)DielectricConstantPT( pressure float64, temperature float64)  (float64, error) {

        p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
        T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
		_, err := if97.GetRegionPT(p, T);
		if err != nil{
			return -1,err
		}
		v,err:=if97.SpecificVolumePT(p, T)
		if err!=nil{
			return -1,err
		}
		retval, err := if97.DielectricConstantRhoT(1 / v, T)
		if err != nil{
			return -1,err
		}
        return retval,nil;

    }

    /**
     * Compression factor (real-gas factor) as a function of pressure &amp;
     * temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return compression factor
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)CompressionFactorPT(pressure float64, temperature float64) (float64,error) {

         p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
         T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature);
         _,err := if97.GetRegionPT(p, T)
		 if err!=nil{
			return -1,err
		 }
		 retval,err:=if97.SpecificVolumePT(p, T)
		 if err!=nil{
			return -1,err
		}

        return 1e3 * p * retval/ (constants.R * T),nil;


    }


    /**
     * Refractive index as a function of pressure, temperature &amp; wave
     * length.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @param wavelength wavelength
     * @return refractive index [-]
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)RefractiveIndexPTLambda(pressure float64, temperature float64,  wavelength float64) (float64,error){

        p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
                T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
                lambda := if97.ConvertToDefault(if97.UnitSystem.WAVELENGTH, wavelength)


			_,err := if97.GetRegionPT(p, T)
			if err!=nil{
				return -1,err
			}
            v,err := if97.SpecificVolumePT(p, T);
			if err!=nil{
				return -1,err
			}

            return if97.RefractiveIndexRhoTLambda(1 / v, T, lambda);

    
    }


    /**
     * Dynamic viscosity as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return dynamic viscosity
     * @throws OutOfRangeException out-of-range exception
     * @see #dynamicViscosityRhoT(double, double)
     */
	 func (if97 *IF97)DynamicViscosityPT(pressure float64, temperature float64)(float64,error) {

        p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
        T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
        var eta float64

 
        _,err := if97.GetRegionPT(p, T)
		if err!=nil{
			return -1,err
		}
		v,err:=if97.SpecificVolumePT(p, T);
		if err!=nil{
			return -1,err
		}
        eta,_ = if97.DynamicViscosityRhoT(1 / v, T);

        return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, eta),nil;
    }

    /**
     * Isothermal compressibility as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return isothermal compressibility
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)CompressibilityPT(pressure float64, temperature float64) (float64,error){

        		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
                T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
                var kappaT float64

 
       _,err := if97.GetRegionPT(p, T)
	   if err!=nil{
		return -1,err
	}
	   
	   kappaT = if97.IsothermalCompressibilityPT(p, T);

        return if97.ConvertFromDefault(if97.UnitSystem.COMPRESSIBILITY, kappaT),nil;
    }


	    /**
     * Heat capacity ratio as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return heat capacity ratio
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)HeatCapacityRatioPT(pressure float64, temperature float64) (float64,error) {

        		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
                T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

				_,err := if97.GetRegionPT(p, T)
				if err!=nil{
				 return -1,err
			 }
		retval:=if97.Region.HeatCapacityRatioPT(p, T);	
			
			
		return if97.ConvertFromDefault(if97.UnitSystem.COMPRESSIBILITY, retval),nil;
    }
	    /**
     * Specific enthalpy as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return specific enthalpy
     * @throws OutOfRangeException out-of-range exception
     */
	func (if97 *IF97) SpecificEnthalpyPT(pressure float64, temperature float64) (float64,error) {

        p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
                T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

				_,err := if97.GetRegionPT(p, T)
				if err!=nil{
				 return -1,err
			 }
  
            h := if97.Region.SpecificEnthalpyPT(p, T);


        return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h),nil;
    }

    /**
     * Isentropic exponent as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return isentropic exponent
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)IsentropicExponentPT(pressure float64, temperature float64) (float64,error) {

		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
		T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

		_,err := if97.GetRegionPT(p, T)
		if err!=nil{
		 	return -1,err
	 	}
		retval:=if97.Region.IsentropicExponentPT(p, T);	

        return  retval,nil;
	 }	
    /**
     * Specific volume as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return specific volume
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)SpecificVolumePT(pressure float64, temperature float64) (float64,error){
		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
		T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

        var v float64

        _,err := if97.GetRegionPT(p, T)
		if err!=nil{
		 	return -1,err
	 	}
			
		v=if97.Region.SpecificVolumePT(p, T);

        return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v),nil;
    }
	    /**
     * Specific entropy as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return specific entropy
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)SpecificEntropyPT(pressure float64, temperature float64) (float64,error){

		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
		T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
    	var s float64
		_,err := if97.GetRegionPT(p, T)
		if err!=nil{
		 	return -1,err
	 	}
   
        s = if97.Region.SpecificEntropyPT(p, T);

    
        return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s),nil;
    }
    func (if97 *IF97)ThermalConductivityRhoT(rho float64, T float64) float64 {

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
	    /**
     * Temperature. [IF97 Supplementary Release S04]
     *
     * @param enthalpy specific enthalpy
     * @param entropy specific entropy
     * @return temperature
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97) TemperatureHS(enthalpy float64, entropy float64) (float64,error) {

		h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
		s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
         var T float64;
		 _,err := if97.GetRegionHS(enthalpy, entropy)
		if err!=nil{
		 	return -1,err
	 	}
      
        T = if97.Region.TemperatureHS(h, s);

    
        return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T),nil;
    }
	
    /**
     * Pressure as a function of specific enthalpy &amp; specific entropy.
     *
     * @param enthalpy specific enthalpy
     * @param entropy specific entropy
     * @return pressure
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97) PressureHS(enthalpy float64, entropy float64) (float64,error) {

        h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
		s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
        var p float64;

        _,err := if97.GetRegionHS(enthalpy, entropy)
		if err!=nil{
		 	return -1,err
	 	}
			
		p = if97.Region.PressureHS(h, s);

    
        return if97.ConvertFromDefault(if97.UnitSystem.PRESSURE, p),nil;
    }

    /**
     * Temperature.
     *
     * @param pressure absolute pressure
     * @param entropy specific entropy
     * @return temperature
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97) TemperaturePS(pressure float64, entropy float64) (float64,error) {

       p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
       s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
    	var T float64;

		_,err := if97.GetRegionPS(p, s)
		if err!=nil{
		 	return -1,err
	 	}
		T=if97.Region.TemperaturePS(p, s)
        return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T),nil;
    }
    /**
     * Temperature.
     *
     * @param pressure absolute pressure [MPa]
     * @param enthalpy specific enthalpy [kJ/(kg)]
     * @return temperature [K]
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97) TemperaturePH(pressure float64, enthalpy float64) (float64,error)  {

        p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
        h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
        var T float64;
		_,err := if97.GetRegionPH(p, h)
		if err!=nil{
			return -1,err
		}
        T=if97.Region.TemperaturePH(p, h)

        return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T),nil;
    }
	   /**
     * Specific volume as a function of pressure &amp; specific enthalpy.
     *
     * @param pressure absolute pressure
     * @param enthalpy specific enthalpy
     * @return specific volume
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97) SpecificVolumePH(pressure float64, enthalpy float64) (float64,error) {

		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
        h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
        var v float64;
		_,err := if97.GetRegionPH(p, h)
		if err!=nil{
			return -1,err
		}
        v = if97.Region.SpecificVolumePH(p, h);

        return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v),nil;
    }
	    /**
     * Specific volume as a function of pressure &amp; specific entropy.
     *
     * @param pressure absolute pressure
     * @param entropy specific entropy
     * @return specific volume
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97) SpecificVolumePS(pressure float64, entropy float64) (float64,error) {


		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
		s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
		 var v float64;
 
		 _,err := if97.GetRegionPS(p, s)
		 if err!=nil{
			  return -1,err
		  }
		 v=if97.Region.SpecificVolumePS(p, s)

        return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v),nil;
    }
	    /**
     * Speed of sound as a function of pressure &amp; specific enthalpy.
     *
     * @param pressure absolute pressure
     * @param enthalpy specific enthalpy
     * @return speed of sound
     * @throws OutOfRangeException out-of-range exception
     * @see #speedOfSoundPT(double, double)
     */
	 func (if97 *IF97) SpeedOfSoundPH(pressure float64, enthalpy float64) (float64,error){
		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
        h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
        var w float64;
		_,err := if97.GetRegionPH(p, h)
		if err!=nil{
			return -1,err
		}
        if (if97.Region.GetName()=="Region 4") {
             w = fourthRegion.REGION4.SpeedOfSoundPH(p, h);
        } else {
            T := if97.Region.TemperaturePH(p, h);
            w = if97.Region.SpeedOfSoundPT(p, T);
		}   
        return if97.ConvertFromDefault(if97.UnitSystem.SPEED_OF_SOUND, w),nil;
    }

	/**
     * Speed of sound as a function of pressure &amp; temperature.
     *
     * @param pressure absolute pressure
     * @param temperature temperature
     * @return speed of sound
     * @throws OutOfRangeException out-of-range exception
     */
	 func (if97 *IF97)SpeedOfSoundPT(pressure float64, temperature float64) (float64,error) {
		p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
		T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
        var w float64;
		_,err := if97.GetRegionPT(p, T)
		if err!=nil{
			return -1,err
		}

       
        w = if97.Region.SpeedOfSoundPT(p, T);

 
        return if97.ConvertFromDefault(if97.UnitSystem.SPEED_OF_SOUND, w),nil;
    }


    func  (if97 *IF97)ThermalDiffusivityPH(p float64, h float64) (float64,error) {
         		v,err  := if97.SpecificVolumePH(p, h)
				if err!=nil{
					return -1,err
				}
				rho:=1/v
                T,err := if97.TemperaturePH(p, h)
				if err!=nil{
					return -1,err
				}
                lambda := if97.ThermalConductivityRhoT(rho, T)
                var cp float64;

        if (if97.Region.GetName()== "Region4") {
            cp = fourthRegion.REGION4.SpecificIsobaricHeatCapacityPH(p, h);

        } else {
            cp = if97.SpecificIsobaricHeatCapacityPT(p, T);
        }
        return lambda / rho / cp,nil
    }