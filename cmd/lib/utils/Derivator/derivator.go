package derivator

import (
	"fmt"

	"if97.com/cmd/lib/fourthRegion"
	"if97.com/cmd/lib/utils/calculator"
	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
)


var (
	ps13 = fourthRegion.SaturationPressureT(constants.T13)
)
/**
 * Gets the partial derivative of z with respect to p_or_T (pressure or
 * temperature) for constant y in SI units.
 *
 * @param pressure pressure [MPa]
 * @param enthalpy specific enthalpy [kJ/kg]
 * @param p_or_T pressure or temperature [IF97.Quantity.p|IF97.Quantity.T]
 * @param y specific enthalpy or specific volume
 * @param z specific enthalpy or specific volume
 * @return partial derivative [SI]
 */
func partialDerivativePH(pressure float64, enthalpy float64, p_or_T quantity.Quantity, y quantity.Quantity, z quantity.Quantity) float64 {

	switch z.NAME {
	case "density":
		rho := 1 / fourthRegion.SpecificVolumePH(pressure, enthalpy) // [kg/m³]

		return -rho * rho * partialDerivativePH(pressure, enthalpy, p_or_T, y, quantity.V)
	}
	var dy, dz float64
	h := specificEnthalpiesP(pressure)          // [kJ/kg]
	v := specificVolumesP(pressure)             // [m³/kg]
	dz_dpT := derivativePQ(pressure, p_or_T, z) // [SI]
	dy_dpT := derivativePQ(pressure, p_or_T, y) // [SI]

	switch z.NAME {
	case "specific enthalpy":
		dz = (h[1] - h[0]) * 1e3 // [J/kg]
		break

	case "specific volume":
		dz = v[1] - v[0] // [m³/kg]
		break

	default:
		panic(fmt.Sprintf("Partial derivative of %s is currently not supported in the saturated region, only specific enthalpy and specific volume.", z))
		//IllegalArgumentException("Partial derivative of " + z + " is currently not supported in the saturated region, only specific enthalpy and specific volume.");
	}
	switch y.NAME {
	case "specific enthalpy":
		dy = (h[0] - h[1]) * 1e3 // [J/kg]
		break

	case "specific volume":
		dy = v[0] - v[1] // [m³/kg]
		break

	default:
		panic(fmt.Sprintf("Partial derivative of %s is currently not supported in the saturated region, only specific enthalpy and specific volume.", y))
		//IllegalArgumentException("Partial derivative for constant " + y + " is currently not supported in the saturated region, only specific enthalpy and specific volume.");
	}
	x := vapourFraction(enthalpy, h) // [-]
	dx_dpT_y := (x*dy_dpT[1] + (1-x)*dy_dpT[0]) / dy

	return dz_dpT[0] + dz*dx_dpT_y + x*(dz_dpT[1]-dz_dpT[0])
}


    /**
     * Gets the derivative of z with respect to pT (pressure or temperature) in
     * SI units, along the saturation line.
     *
     * @param pressure [MPa]
     * @param pT pressure or temperature [IF97.Quantity.p | IF97.Quantity.T]
     * @param z any quantity
     * @return derivative (dz/dp)T or (dz/dT)p [SI]
     */
  func derivativePQ(pressure float64,  pT quantity.Quantity,  z  quantity.Quantity) []float64 {

         T := fourthRegion.SaturationTemperatureP(pressure) // [K]
                dT_dp := derivativeP(pressure) // [K/Pa]
         var dz_dp_T [2]float64
                var dz_dT_p  [2]float64

        if ps13 <= pressure {
            /*
             Region 3
             */
             h := specificEnthalpiesP(pressure) // [kJ/kg]
                    rho := densitiesRegion3(pressure, h) // [kg/m³]
            dz_dp_T[0] = calculator.PartialDerivativeRhoT(rho[0], T, quantity.P, quantity.T, z); // [SI]
            dz_dp_T[1] = calculator.PartialDerivativeRhoT(rho[1], T, quantity.P, quantity.T, z); // [SI]
            dz_dT_p[0] = calculator.PartialDerivativeRhoT(rho[0], T, quantity.T, quantity.P, z); // [SI]
            dz_dT_p[1] = calculator.PartialDerivativeRhoT(rho[1], T, quantity.T, quantity.P, z); // [SI]

        } else {
            /*
             Regions 1 & 2
             */
            dz_dp_T[0] = calculator.PartialDerivativePT(Region.REGION1, pressure, T, quantity.P, quantity.T, z); // [SI]
            dz_dp_T[1] = calculator.PartialDerivativePT(Region.REGION2, pressure, T, quantity.P, quantity.T, z); // [SI]
            dz_dT_p[0] = calculator.PartialDerivativePT(Region.REGION1, pressure, T, quantity.T, quantity.P, z); // [SI]
            dz_dT_p[1] = calculator.PartialDerivativePT(Region.REGION2, pressure, T, quantity.T, quantity.P, z); // [SI]
        }
        switch (pT.SYMBOL) {
            case "p":
                return []float64{
                    dz_dp_T[0] + dz_dT_p[0] * dT_dp,
                    dz_dp_T[1] + dz_dT_p[1] * dT_dp};

            case "T":
                return []float64{
                    dz_dT_p[0] + dz_dp_T[0] / dT_dp,
                    dz_dT_p[1] + dz_dp_T[1] / dT_dp};

            default:
                panic(fmt.Sprintf("Quantity pT should be either pressure or temperature, not:%s " , pT)) //IllegalArgumentException("Quantity pT should be either pressure or temperature, not: " + pT);
        }
    }