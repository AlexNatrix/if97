package IF97

import (
	"errors"
	"math"

	"if97.com/cmd/lib/fourthRegion"
	rangeError "if97.com/cmd/lib/region/range_error"
	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
)

var p0 = fourthRegion.REGION4.SaturationPressureT(constants.T0)

/**
 * Isothermal compressibility as a function of pressure &amp; specific
 * entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return isothermal compressibility
 * @return RangeError
 */
func (if97 *IF97) CompressibilityPS(pressure float64, entropy float64) (float64, error) {
	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var kappaT float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	kappaT = if97.region.IsothermalCompressibilityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.COMPRESSIBILITY, kappaT), nil
}




/**
 * Isothermal compressibility as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return isothermal compressibility
 * @return RangeError
 */
func (if97 *IF97) CompressibilityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var kappaT float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}

	kappaT = if97.region.IsothermalCompressibilityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.COMPRESSIBILITY, kappaT), nil
}

/**
 * Isothermal compressibility as a function of pressure &amp; specific
 * enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return isothermal compressibility
 * @return RangeError
 */
func (if97 *IF97) CompressibilityPH(pressure float64, enthalpy float64) (float64, error) {
	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var kappaT float64
	reg, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		//TODO bad code
		kappaT = fourthRegion.REGION4.IsothermalCompressibilityPH(p, h)

	} else {
		T := if97.region.TemperaturePH(p, h)

		kappaT = if97.region.IsothermalCompressibilityPT(p, T)
	}
	return if97.ConvertFromDefault(if97.UnitSystem.COMPRESSIBILITY, kappaT), nil
}

/**
 * Isothermal compressibility as a function of specific enthalpy &amp;
 * specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return isothermal compressibility
 * @return RangeError
 */
func (if97 *IF97) CompressibilityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var kappaT float64
	_, err := if97.region.GetRegionHS(enthalpy, entropy)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	kappaT = if97.region.IsothermalCompressibilityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.COMPRESSIBILITY, kappaT), nil
}

/**
 * Density as a function of specific enthalpy &amp; specific entropy.
 *
 * <p>
 * This is a convenience method which simply calls
 * <code>1.0 / specificVolumeHS(enthalpy, entropy)</code>.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return density
 * @return RangeError
 * @see #specificVolumeHS(double, double)
 */
func (if97 *IF97) DensityHS(enthalpy float64, entropy float64) (float64, error) {
	ret, err := if97.SpecificVolumeHS(enthalpy, entropy)

	if err != nil {
		return -1, err
	}

	return 1.0 / ret, err
}

/**
 * Density as a function of pressure &amp; specific enthalpy.
 *
 * <p>
 * This is a convenience method which simply calls
 * <code>1.0 / specificVolumePH(pressure, enthalpy)</code>.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return density
 * @return RangeError
 * @see #specificVolumePH(double, double)
 */
func (if97 *IF97) DensityPH(pressure float64, enthalpy float64) (float64, error) {
	ret, err := if97.SpecificVolumePH(pressure, enthalpy)

	if err != nil {
		return -1, err
	}

	return 1.0 / ret, err
}

/**
 * Density as a function of pressure &amp; specific entropy.
 *
 * <p>
 * This is a convenience method which simply calls
 * <code>1.0 / specificVolumePS(pressure, entropy)</code>.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return density
 * @return RangeError
 * @see #specificVolumePS(double, double)
 */
func (if97 *IF97) DensityPS(pressure float64, entropy float64) (float64, error) {
	ret, err := if97.SpecificVolumePS(pressure, entropy)
	if err != nil {
		return -1, err
	}

	return 1.0 / ret, err
}

/**
 * Density as a function of pressure &amp; temperature.
 *
 * <p>
 * This is a convenience method which simply calls
 * <code>1.0 / specificVolumePT(pressure, temperature)</code>.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return density
 * @return RangeError
 * @see #specificVolumePT(double, double)
 */
func (if97 *IF97) DensityPT(pressure float64, temperature float64) (float64, error) {
	ret, err := if97.SpecificVolumePT(pressure, temperature)

	if err != nil {
		return -1, err
	}

	return 1.0 / ret, err
}

/**
 * Density as a function of pressure &amp; vapour fraction.
 *
 * <p>
 * This is a convenience method which simply calls
 * <code>1.0 / specificVolumePX(pressure, vapour fraction)</code>.
 *
 * @param pressure absolute pressure
 * @param vapourFraction vapour fraction
 * @return density
 * @return RangeError
 * @see #specificVolumePX(double, double)
 */
func (if97 *IF97) DensityPX(pressure float64, vapourFraction float64) (float64, error) {
	ret, err := if97.SpecificVolumePX(pressure, vapourFraction)

	if err != nil {
		return -1, err
	}

	return 1.0 / ret, err
}

/**
 * Density as a function of temperature &amp; vapour fraction.
 *
 * <p>
 * This is a convenience method which simply calls
 * <code>1.0 / specificVolumeTX(temperature, vapour fraction)</code>.
 *
 * @param temperature temperature
 * @param vapourFraction vapour fraction
 * @return density
 * @return RangeError
 * @see #specificVolumeTX(double, double)
 */
func (if97 *IF97) DensityTX(temperature float64, vapourFraction float64) (float64, error) {
	ret, err := if97.SpecificVolumeTX(temperature, vapourFraction)

	if err != nil {
		return -1, err
	}

	return 1.0 / ret, err
}

/**
 * Dielectric constant as a function of specific enthalpy and specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return dielectric constant [-]
 * @return RangeError
 * @see #dielectricConstantRhoT(double, double)
 */
func (if97 *IF97) DielectricConstantHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumeHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	return if97.DielectricConstantRhoT(1.0/v, T)

}

/**
 * Dielectric constant as a function of pressure and specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return dielectric constant [-]
 * @return RangeError
 * @see #dielectricConstantRhoT(double, double)
 */
func (if97 *IF97) DielectricConstantPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePH(p, h)
	T := if97.region.TemperaturePH(p, h)

	return if97.DielectricConstantRhoT(1.0/v, T)
}

/**
 * Dielectric constant as a function of pressure and specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return dielectric constant [-]
 * @return RangeError
 * @see #dielectricConstantRhoT(double, double)
 */
func (if97 *IF97) DielectricConstantPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	_, err := if97.region.GetRegionPS(p, s)

	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumePS(p, s)
	T := if97.region.TemperaturePS(p, s)

	return if97.DielectricConstantRhoT(1.0/v, T)
}

/**
 * Dielectric constant as a function of pressure and temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return dielectric constant [-]
 * @return RangeError
 * @see #dielectricConstantRhoT(double, double)
 */
func (if97 *IF97) DielectricConstantPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePT(p, T)

	return if97.DielectricConstantRhoT(1.0/v, T)

}

/**
 * Dielectric constant (relative static dielectric constant or relative
 * static permittivity) as a function of density and temperature.
 *
 * @param density density
 * @param temperature temperature
 * @return dielectric constant [-]
 * @return RangeError
 */
func (if97 *IF97) DielectricConstantRhoT(density float64, temperature float64) (float64, error) {

	rho := if97.ConvertToDefault(if97.UnitSystem.DENSITY, density)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

	return if97.defaultDielectricConstantRhoT(rho, T)

}

/**
 * Dynamic viscosity as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return dynamic viscosity
 * @return RangeError
 * @see #dynamicViscosityRhoT(double, double)
 */
func (if97 *IF97) DynamicViscosityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var eta float64
	_, err := if97.region.GetRegionHS(h, s)

	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumeHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	eta, err = if97.DynamicViscosityRhoT(1.0/v, T)

	if err != nil {
		return -1, err
	}

	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, eta), err
}

/**
 * Dynamic viscosity as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return dynamic viscosity
 * @return RangeError
 * @see #dynamicViscosityRhoT(double, double)
 */
func (if97 *IF97) DynamicViscosityPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var eta float64
	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePH(p, h)
	T := if97.region.TemperaturePH(p, h)

	eta, err = if97.DynamicViscosityRhoT(1.0/v, T)

	if err != nil {
		return -1, err
	}

	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, eta), err
}

/**
 * Dynamic viscosity as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return dynamic viscosity
 * @return RangeError
 * @see #dynamicViscosityRhoT(double, double)
 */
func (if97 *IF97) DynamicViscosityPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var eta float64
	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumePS(p, s)
	T := if97.region.TemperaturePS(p, s)

	eta, err = if97.DynamicViscosityRhoT(1.0/v, T)
	if err != nil {
		return -1, err
	}
	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, eta), err
}

/**
 * Dynamic viscosity as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return dynamic viscosity
 * @return RangeError
 * @see #dynamicViscosityRhoT(double, double)
 */
func (if97 *IF97) DynamicViscosityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var eta float64
	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumePT(p, T)

	eta, err = if97.DynamicViscosityRhoT(1.0/v, T)
	if err != nil {
		return -1, err
	}
	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, eta), err
}

/**
 * Dynamic viscosity as a function of density &amp; temperature.
 *
 * @param density density
 * @param temperature temperature
 * @return dynamic viscosity
 * @return RangeError
 */
func (if97 *IF97) DynamicViscosityRhoT(density float64, temperature float64) (float64, error) {

	rho := if97.ConvertToDefault(if97.UnitSystem.DENSITY, density)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	eta, err := if97.defaultDynamicViscosityRhoT(rho, T)
	if err != nil {
		return -1, err
	}
	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, eta), err
}

// UTILITY
func (if97 *IF97) getRegionPT(pressure float64, temperature float64) (int, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

	return if97.region.GetRegionPT(p, T)
}

/**
 * Heat capacity ratio as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return heat capacity ratio
 * @return RangeError
 * @see #heatCapacityRatioPT(double, double)
 */
func (if97 *IF97) HeatCapacityRatioHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	return if97.region.HeatCapacityRatioPT(p, T), err
}

/**
 * Heat capacity ratio as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return heat capacity ratio
 * @return RangeError
 * @see #heatCapacityRatioPT(double, double)
 */
func (if97 *IF97) HeatCapacityRatioPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)

	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePH(p, h)

	return if97.region.HeatCapacityRatioPT(p, T), err
}

/**
 * Heat capacity ratio as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return heat capacity ratio
 * @return RangeError
 * @see #heatCapacityRatioPT(double, double)
 */
func (if97 *IF97) HeatCapacityRatioPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	return if97.region.HeatCapacityRatioPT(p, T), err
}

/**
 * Heat capacity ratio as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return heat capacity ratio
 * @return RangeError
 */
func (if97 *IF97) HeatCapacityRatioPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}

	return if97.region.HeatCapacityRatioPT(p, T), err

}

/**
 * Isentropic exponent as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return isentropic exponent
 * @return RangeError
 * @see #isentropicExponentPT(double, double)
 */
func (if97 *IF97) IsentropicExponentHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	return if97.region.IsentropicExponentPT(p, T), err

}

/**
 * Isentropic exponent as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return isentropic exponent
 * @return RangeError
 * @see #isentropicExponentPT(double, double)
 */
func (if97 *IF97) IsentropicExponentPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePH(p, h)

	return if97.region.IsentropicExponentPT(p, T), err

}

/**
 * Isentropic exponent as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return isentropic exponent
 * @return RangeError
 * @see #isentropicExponentPT(double, double)
 */
func (if97 *IF97) IsentropicExponentPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	return if97.region.IsentropicExponentPT(p, T), err

}

/**
 * Isentropic exponent as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return isentropic exponent
 * @return RangeError
 */
func (if97 *IF97) IsentropicExponentPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}

	return if97.region.IsentropicExponentPT(p, T), err

}

/**
 * Isobaric cubic expansion coefficient as a function of specific enthalpy
 * &amp; specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return isobaric cubic expansion coefficient
 * @return RangeError
 * @see #isobaricCubicExpansionCoefficientPT(double, double)
 */
func (if97 *IF97) IsobaricCubicExpansionCoefficientHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var alphaV float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	alphaV = if97.region.IsobaricCubicExpansionCoefficientPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.ISOBARIC_CUBIC_EXPANSION_COEFFICIENT, alphaV), err
}

/**
 * Isobaric cubic expansion coefficient as a function of pressure &amp;
 * specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return isobaric cubic expansion coefficient
 * @return RangeError
 * @see #isobaricCubicExpansionCoefficientPT(double, double)
 */
func (if97 *IF97) IsobaricCubicExpansionCoefficientPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var alphaV float64

	reg, err := if97.region.GetRegionPS(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		alphaV = fourthRegion.REGION4.IsobaricCubicExpansionCoefficientPH(p, h)

	} else {
		T := if97.region.TemperaturePH(p, h)

		alphaV = if97.region.IsobaricCubicExpansionCoefficientPT(p, T)
	}

	return if97.ConvertFromDefault(if97.UnitSystem.ISOBARIC_CUBIC_EXPANSION_COEFFICIENT, alphaV), err
}

/**
 * Isobaric cubic expansion coefficient as a function of pressure &amp;
 * specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return isobaric cubic expansion coefficient
 * @return RangeError
 * @see #isobaricCubicExpansionCoefficientPT(double, double)
 */
func (if97 *IF97) IsobaricCubicExpansionCoefficientPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var alphaV float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	alphaV = if97.region.IsobaricCubicExpansionCoefficientPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.ISOBARIC_CUBIC_EXPANSION_COEFFICIENT, alphaV), err
}

/**
 * Isobaric cubic expansion coefficient as a function of pressure &amp;
 * temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return isobaric cubic expansion coefficient
 * @return RangeError
 */
func (if97 *IF97) IsobaricCubicExpansionCoefficientPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var alphaV float64
	_, err := if97.region.GetRegionPS(p, T)
	if err != nil {
		return -1, err
	}

	alphaV = if97.region.IsobaricCubicExpansionCoefficientPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.ISOBARIC_CUBIC_EXPANSION_COEFFICIENT, alphaV), err
}

/**
 * Isobaric cubic expansion coefficient as a function of pressure &amp;
 * vapour fraction.
 *
 * @param pressure absolute pressure
 * @param vapourFraction vapour fraction
 * @return isobaric cubic expansion coefficient
 * @return RangeError
 * @see #isobaricCubicExpansionCoefficientPH(double, double)
 */
func (if97 *IF97) IsobaricCubicExpansionCoefficientPX(pressure float64, vapourFraction float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	h := fourthRegion.REGION4.SpecificEnthalpyPX(p, vapourFraction)
	alphaV := fourthRegion.REGION4.IsobaricCubicExpansionCoefficientPH(p, h)

	return if97.ConvertFromDefault(if97.UnitSystem.ISOBARIC_CUBIC_EXPANSION_COEFFICIENT, alphaV), err
}

/**
 * Isobaric cubic expansion coefficient as a function of temperature &amp;
 * vapour fraction.
 *
 * @param temperature temperature
 * @param vapourFraction vapour fraction
 * @return isobaric cubic expansion coefficient
 * @return RangeError
 * @see #isobaricCubicExpansionCoefficientPH(double, double)
 */
func (if97 *IF97) IsobaricCubicExpansionCoefficientTX(temperature float64, vapourFraction float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	p := fourthRegion.REGION4.SaturationPressureT(T)
	h := fourthRegion.REGION4.SpecificEnthalpyPX(p, vapourFraction)
	alphaV := fourthRegion.REGION4.IsobaricCubicExpansionCoefficientPH(p, h)

	return if97.ConvertFromDefault(if97.UnitSystem.ISOBARIC_CUBIC_EXPANSION_COEFFICIENT, alphaV), err
}

/**
 * Specific isobaric heat capacity as a function of specific enthalpy &amp;
 * specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return specific isobaric heat capacity
 * @return RangeError
 * @see #isobaricHeatCapacityPT(double, double)
 */
func (if97 *IF97) IsobaricHeatCapacityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var cp float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	cp = if97.region.SpecificIsobaricHeatCapacityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cp), err
}

/**
* Specific isobaric heat capacity as a function of pressure &amp; specific
* enthalpy.
*
* @param pressure absolute pressure
* @param enthalpy specific enthalpy
* @return specific isobaric heat capacity
* @return RangeError
 * @see #isobaricHeatCapacityPT(double, double)
*/
func (if97 *IF97) IsobaricHeatCapacityPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var cp float64

	reg, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		cp = fourthRegion.REGION4.SpecificIsobaricHeatCapacityPH(p, h)

	} else {
		T := if97.region.TemperaturePH(p, h)

		cp = if97.region.SpecificIsobaricHeatCapacityPT(p, T)
	}

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cp), err
}

/**
 * Specific isobaric heat capacity as a function of pressure &amp; specific
 * entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return specific isobaric heat capacity
 * @return RangeError
 * @see #isobaricHeatCapacityPT(double, double)
 */
func (if97 *IF97) IsobaricHeatCapacityPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var cp float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	cp = if97.region.SpecificIsobaricHeatCapacityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cp), err
}

/**
 * Specific isobaric heat capacity as a function of pressure &amp;
 * temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific isobaric heat capacity
 * @return RangeError
 */
func (if97 *IF97) IsobaricHeatCapacityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var cp float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	cp = if97.region.SpecificIsobaricHeatCapacityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cp), err
}

/**
 * Specific isochoric heat capacity as a function of specific enthalpy &amp;
 * specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return specific isochoric heat capacity
 * @return RangeError
 * @see #isochoricHeatCapacityPT(double, double)
 */
func (if97 *IF97) IsochoricHeatCapacityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var cv float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	cv = if97.region.SpecificIsochoricHeatCapacityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cv), err
}

/**
 * Specific isochoric heat capacity as a function of pressure &amp; specific
 * enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return specific isochoric heat capacity
 * @return RangeError
 * @see #isochoricHeatCapacityPT(double, double)
 */
func (if97 *IF97) IsochoricHeatCapacityPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var cv float64

	reg, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		cv = fourthRegion.REGION4.SpecificIsochoricHeatCapacityPH(p, h)

	} else {
		T := if97.region.TemperaturePH(p, h)

		cv = if97.region.SpecificIsochoricHeatCapacityPT(p, T)
	}

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cv), err
}

/**
 * Specific isochoric heat capacity as a function of pressure &amp; specific
 * entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return specific isochoric heat capacity
 * @return RangeError
 * @see #isochoricHeatCapacityPT(double, double)
 */
func (if97 *IF97) IsochoricHeatCapacityPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var cv float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	cv = if97.region.SpecificIsochoricHeatCapacityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cv), err
}

/**
 * Specific isochoric heat capacity as a function of pressure &amp;
 * temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific isochoric heat capacity
 * @return RangeError
 */
func (if97 *IF97) IsochoricHeatCapacityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var cv float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}

	cv = if97.region.SpecificIsochoricHeatCapacityPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_HEAT_CAPACITY, cv), err
}

/**
 * Kinematic viscosity as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return kinematic viscosity
 * @return RangeError
 * @see #kinematicViscosityRhoT(double, double)
 */
func (if97 *IF97) KinematicViscosityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var nu float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumeHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	eps, err := if97.defaultDynamicViscosityRhoT(1.0/v, T)
	if err != nil {
		return -1, err
	}
	nu = eps * v
	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, nu), err
}

/**
 * Kinematic viscosity as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return kinematic viscosity
 * @return RangeError
 * @see #kinematicViscosityRhoT(double, double)
 */
func (if97 *IF97) KinematicViscosityPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var nu float64

	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePH(p, h)
	T := if97.region.TemperaturePH(p, h)

	eta, err := if97.defaultDynamicViscosityRhoT(1.0/v, T)
	if err != nil {
		return -1, err
	}
	nu = eta * v

	return if97.ConvertFromDefault(if97.UnitSystem.KINEMATIC_VISCOSITY, nu), err
}

/**
 * Kinematic viscosity as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return kinematic viscosity
 * @return RangeError
 * @see #kinematicViscosityRhoT(double, double)
 */
func (if97 *IF97) KinematicViscosityPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var nu float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePS(p, s)
	T := if97.region.TemperaturePS(p, s)

	eta, err := if97.defaultDynamicViscosityRhoT(1.0/v, T)
	if err != nil {
		return -1, err
	}
	nu = eta * v

	return if97.ConvertFromDefault(if97.UnitSystem.DYNAMIC_VISCOSITY, nu), err
}

/**
 * Kinematic viscosity as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return kinematic viscosity
 * @return RangeError
 * @see #kinematicViscosityRhoT(double, double)
 */
func (if97 *IF97) KinematicViscosityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var nu float64

	_, err := if97.region.GetRegionPS(p, T)
	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumePT(p, T)

	eta, err := if97.defaultDynamicViscosityRhoT(1.0/v, T)
	if err != nil {
		return -1, err
	}
	nu = eta * v

	return if97.ConvertFromDefault(if97.UnitSystem.KINEMATIC_VISCOSITY, nu), err
}

/**
 * Kinematic viscosity as a function of density &amp; temperature.
 *
 * @param density density
 * @param temperature temperature
 * @return kinematic viscosity
 * @return RangeError
 */
func (if97 *IF97) KinematicViscosityRhoT(density float64, temperature float64) (float64, error) {

	rho := if97.ConvertToDefault(if97.UnitSystem.DENSITY, density)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var nu float64

	eta, err := if97.defaultDynamicViscosityRhoT(rho, T)
	if err != nil {
		return -1, err
	}
	nu = eta / rho

	return if97.ConvertFromDefault(if97.UnitSystem.KINEMATIC_VISCOSITY, nu), err
}

/**
 * Pressure as a function of specific enthalpy &amp; specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return pressure
 * @return RangeError
 */
func (if97 *IF97) PressureHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var p float64

	reg, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}
	p = if97.region.PressureHS(h, s)
	if reg==5 && math.IsNaN(p){
		return -1, errors.New("region 5 PressureHS is unimplemented, Im sorry")
	}
	return if97.ConvertFromDefault(if97.UnitSystem.PRESSURE, p), err
}

/**
 * Refractive index as a function of specific enthalpy, specific entropy
 * &amp; wave length.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @param wavelength wavelength
 * @return refractive index [-]
 * @return RangeError
 */
func (if97 *IF97) RefractiveIndexHSLambda(enthalpy float64, entropy float64, wavelength float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	lambda := if97.ConvertToDefault(if97.UnitSystem.WAVELENGTH, wavelength)

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumeHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	return defaultRefractiveIndexRhoTLambda(1.0/v, T, lambda)

}

/**
 * Refractive index as a function of pressure, specific enthalpy &amp; wave
 * length.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @param wavelength wavelength
 * @return refractive index [-]
 * @return RangeError
 */
func (if97 *IF97) RefractiveIndexPHLambda(pressure float64, enthalpy float64, wavelength float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	lambda := if97.ConvertToDefault(if97.UnitSystem.WAVELENGTH, wavelength)

	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePH(p, h)
	T := if97.region.TemperaturePH(p, h)

	return defaultRefractiveIndexRhoTLambda(1/v, T, lambda)

}

/**
 * Refractive index as a function of pressure, specific entropy &amp; wave
 * length.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @param wavelength wavelength
 * @return refractive index [-]
 * @return RangeError
 */
func (if97 *IF97) RefractiveIndexPSLambda(pressure float64, entropy float64, wavelength float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	lambda := if97.ConvertToDefault(if97.UnitSystem.WAVELENGTH, wavelength)

	if _, err := if97.region.GetRegionPS(p, s); err != nil {
		return -1, err
	}

	v := if97.region.SpecificVolumePS(p, s)
	T := if97.region.TemperaturePS(p, s)

	return defaultRefractiveIndexRhoTLambda(1/v, T, lambda)
}

/**
 * Refractive index as a function of pressure, temperature &amp; wave
 * length.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @param wavelength wavelength
 * @return refractive index [-]
 * @return RangeError
 */
func (if97 *IF97) RefractiveIndexPTLambda(pressure float64, temperature float64, wavelength float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	lambda := if97.ConvertToDefault(if97.UnitSystem.WAVELENGTH, wavelength)

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	v := if97.region.SpecificVolumePT(p, T)

	return defaultRefractiveIndexRhoTLambda(1/v, T, lambda)

}

/**
 * Refractive index as a function of density, temperature &amp; wave length.
 *
 * @param density density
 * @param temperature temperature
 * @param waveLength wave length
 * @return refractive index [-]
 * @return RangeError
 */
func (if97 *IF97) RefractiveIndexRhoTLambda(density float64, temperature float64, wavelength float64) (float64, error) {

	rho := if97.ConvertToDefault(if97.UnitSystem.DENSITY, density)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	lambda := if97.ConvertToDefault(if97.UnitSystem.WAVELENGTH, wavelength)

	return defaultRefractiveIndexRhoTLambda(rho, T, lambda)
}

// /**
//   - Saturation pressure as a function of specific enthalpy &amp; specific
//   - entropy.
//     *
//   - @param enthalpy specific enthalpy
//   - @param entropy specific entropy
//   - @return saturation pressure
//   - @return RangeError
//     */
func (if97 *IF97) SaturationPressureHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var p float64

	err := fourthRegion.CheckHS(h, s)
	if err != nil {
		return -1, err
	}

	p = fourthRegion.REGION4.PressureHS(h, s)

	if p < p0 {
		return -1, rangeError.ErrorFromValue(quantity.P, p, p0)
	}
	return if97.ConvertFromDefault(if97.UnitSystem.PRESSURE, p), err
}

/**
 * Saturation pressure as a function of temperature.
 *
 * @param temperature saturation temperature
 * @return saturation pressure
 * @return RangeError
 */
func (if97 *IF97) SaturationPressureT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var p float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	p = fourthRegion.REGION4.SaturationPressureT(T)

	return if97.ConvertFromDefault(if97.UnitSystem.PRESSURE, p), err
}

/**
 * Saturation temperature as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return saturation temperature
 * @return RangeError
 */
func (if97 *IF97) SaturationTemperatureHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var T float64

	err := fourthRegion.CheckHS(h, s)
	if err != nil {
		return -1, err
	}

	T = fourthRegion.REGION4.TemperatureHS(h, s)

	return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T), err
}

/**
 * Saturation temperature as a function of pressure.
 *
 * @param pressure saturation pressure
 * @return saturation temperature
 * @return RangeError
 */
func (if97 *IF97) SaturationTemperatureP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var T float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	T = fourthRegion.REGION4.SaturationTemperatureP(p)
	return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T), err
}

/**
 * Specific enthalpy as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return specific enthalpy
 * @return RangeError
 */
func (if97 *IF97) SpecificEnthalpyPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var h float64
	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	h = if97.region.SpecificEnthalpyPS(p, s)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific enthalpy
 * @return RangeError
 */
func (if97 *IF97) SpecificEnthalpyPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var h float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	h = if97.region.SpecificEnthalpyPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of pressure &amp; vapour fraction.
 *
 * @param pressure absolute pressure
 * @param vapourFraction vapour fraction [-]
 * @return specific enthalpy
 * @return RangeError
 */
func (if97 *IF97) SpecificEnthalpyPX(pressure float64, vapourFraction float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var h float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	h = fourthRegion.REGION4.SpecificEnthalpyPX(p, vapourFraction)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of pressure for saturated liquid.
 *
 * @param pressure saturation pressure
 * @return specific enthalpy
 * @return RangeError
 * @see #specificEnthalpyPX(double, double)
 */
func (if97 *IF97) SpecificEnthalpySaturatedLiquidP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var h float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	h = fourthRegion.REGION4.SpecificEnthalpySaturatedLiquidP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of temperature for saturated liquid.
 *
 * @param temperature saturation temperature
 * @return specific enthalpy
 * @return RangeError
 * @see #specificEnthalpyTX(double, double)
 */
func (if97 *IF97) SpecificEnthalpySaturatedLiquidT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var h float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	h = fourthRegion.REGION4.SpecificEnthalpySaturatedLiquidP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of pressure for saturated vapour.
 *
 * @param pressure saturation pressure
 * @return specific enthalpy
 * @return RangeError
 * @see #specificEnthalpyPX(double, double)
 */
func (if97 *IF97) SpecificEnthalpySaturatedVapourP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var h float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	h = fourthRegion.REGION4.SpecificEnthalpySaturatedVapourP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of temperature for saturated vapour.
 *
 * @param temperature saturation temperature
 * @return specific enthalpy
 * @return RangeError
 * @see #specificEnthalpyTX(double, double)
 */
func (if97 *IF97) SpecificEnthalpySaturatedVapourT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var h float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	h = fourthRegion.REGION4.SpecificEnthalpySaturatedVapourP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific enthalpy as a function of temperature &amp; vapour fraction.
 *
 * @param temperature temperature
 * @param vapourFraction vapour fraction [-]
 * @return specific enthalpy
 * @return RangeError
 */
func (if97 *IF97) SpecificEnthalpyTX(temperature float64, vapourFraction float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var h float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	h = fourthRegion.REGION4.SpecificEnthalpyPX(
		fourthRegion.REGION4.SaturationPressureT(T), vapourFraction,
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, h), err
}

/**
 * Specific entropy as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return specific entropy
 * @return RangeError
 */
func (if97 *IF97) SpecificEntropyPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var s float64

	reg, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		s = fourthRegion.REGION4.SpecificEntropyPH(p, h)

	} else {
		T := if97.region.TemperaturePH(p, h)

		s = if97.region.SpecificEntropyPT(p, T)
	}
	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific entropy
 * @return RangeError
 */
func (if97 *IF97) SpecificEntropyPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var s float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	s = if97.region.SpecificEntropyPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of pressure &amp; vapour fraction.
 *
 * @param pressure absolute pressure
 * @param vapourFraction vapour fraction [-]
 * @return specific entropy
 * @return RangeError
 */
func (if97 *IF97) SpecificEntropyPX(pressure float64, vapourFraction float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var s float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificEntropyPX(p, vapourFraction)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of pressure for saturated liquid.
 *
 * @param pressure saturation pressure
 * @return specific entropy
 * @return RangeError
 * @see #specificEntropyPX(double, double)
 */
func (if97 *IF97) SpecificEntropySaturatedLiquidP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var s float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificEntropySaturatedLiquidP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of temperature for saturated liquid.
 *
 * @param temperature saturation temperature
 * @return specific entropy
 * @return RangeError
 * @see #specificEntropyTX(double, double)
 */
func (if97 *IF97) SpecificEntropySaturatedLiquidT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var s float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificEntropySaturatedLiquidP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of pressure for saturated vapour.
 *
 * @param pressure saturation pressure
 * @return specific entropy
 * @return RangeError
 * @see #specificEntropyPX(double, double)
 */
func (if97 *IF97) SpecificEntropySaturatedVapourP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var s float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificEntropySaturatedVapourP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of temperature for saturated vapour.
 *
 * @param temperature saturation temperature
 * @return specific entropy
 * @return RangeError
 * @see #specificEntropyTX(double, double)
 */
func (if97 *IF97) SpecificEntropySaturatedVapourT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var s float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificEntropySaturatedVapourP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific entropy as a function of temperature &amp; vapour fraction.
 *
 * @param temperature temperature
 * @param vapourFraction vapour fraction [-]
 * @return specific entropy
 * @return RangeError
 */
func (if97 *IF97) SpecificEntropyTX(temperature float64, vapourFraction float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var s float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificEntropyPX(
		fourthRegion.REGION4.SaturationPressureT(T), vapourFraction,
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific Gibbs free energy as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific Gibbs free energy
 * @return RangeError
 */
func (if97 *IF97) SpecificGibbsFreeEnergyPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var g float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	g = if97.region.SpecificGibbsFreeEnergyPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENERGY, g), err
}

/**
 * Specific internal energy as a function of specific enthalpy &amp;
 * specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyPT(double, double)
 */
func (if97 *IF97) SpecificInternalEnergyHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var u float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	u = if97.region.SpecificInternalEnergyPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENERGY, u), err
}

/**
 * Specific internal energy as a function of pressure &amp; specific
 * enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyPT(double, double)
 */
func (if97 *IF97) SpecificInternalEnergyPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var u float64

	reg, err := if97.region.GetRegionPS(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		u = fourthRegion.REGION4.SpecificInternalEnergyPH(p, h)

	} else {
		T := if97.region.TemperaturePH(p, h)

		u = if97.region.SpecificInternalEnergyPT(p, T)
	}

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENERGY, u), err
}

/**
 * Specific internal energy as a function of pressure &amp; specific
 * entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyPT(double, double)
 */
func (if97 *IF97) SpecificInternalEnergyPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var u float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	u = if97.region.SpecificInternalEnergyPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENERGY, u), err
}

/**
 * Specific internal energy as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific internal energy
 * @return RangeError
 */
func (if97 *IF97) SpecificInternalEnergyPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var u float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	u = if97.region.SpecificInternalEnergyPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENERGY, u), err
}

/**
 * Specific internal energy as a function of pressure &amp; vapour fraction.
 *
 * @param pressure absolute pressure
 * @param vapourFraction vapour fraction [-]
 * @return specific internal energy
 * @return RangeError
 */
func (if97 *IF97) SpecificInternalEnergyPX(pressure float64, vapourFraction float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var v float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificInternalEnergyPX(p, vapourFraction)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific internal energy as a function of pressure for saturated liquid.
 *
 * @param pressure saturation pressure
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyPX(double, double)
 */
func (if97 *IF97) SpecificInternalEnergySaturatedLiquidP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var s float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificInternalEnergySaturatedLiquidP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific internal energy as a function of temperature for saturated
 * liquid.
 *
 * @param temperature saturation temperature
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyTX(double, double)
 */
func (if97 *IF97) SpecificInternalEnergySaturatedLiquidT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var s float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificInternalEnergySaturatedLiquidP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific internal energy as a function of pressure for saturated vapour.
 *
 * @param pressure saturation pressure
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyPX(double, double)
 */
func (if97 *IF97) SpecificInternalEnergySaturatedVapourP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var s float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificInternalEnergySaturatedVapourP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific internal energy as a function of temperature for saturated
 * vapour.
 *
 * @param temperature saturation temperature
 * @return specific internal energy
 * @return RangeError
 * @see #specificInternalEnergyTX(double, double)
 */
func (if97 *IF97) SpecificInternalEnergySaturatedVapourT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var s float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	s = fourthRegion.REGION4.SpecificInternalEnergySaturatedVapourP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_ENTROPY, s), err
}

/**
 * Specific internal energy as a function of temperature &amp; vapour
 * fraction.
 *
 * @param temperature temperature
 * @param vapourFraction vapour fraction [-]
 * @return specific internal energy
 * @return RangeError
 */
func (if97 *IF97) SpecificInternalEnergyTX(temperature float64, vapourFraction float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var v float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificInternalEnergyPX(
		fourthRegion.REGION4.SaturationPressureT(T), vapourFraction,
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumeHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var v float64
	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	v = if97.region.SpecificVolumeHS(h, s)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumePH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var v float64

	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}
	v = if97.region.SpecificVolumePH(p, h)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumePS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var v float64
	reg, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}
	v = if97.region.SpecificVolumePS(p, s)
	if math.IsNaN(v) && reg ==3 {
		return -1, errors.New("unsupported subregion in region 3 SpecificVolumePS")
	}
	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumePT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var v float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	v = if97.region.SpecificVolumePT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of pressure &amp; vapour fraction.
 *
 * @param pressure absolute pressure
 * @param vapourFraction vapour fraction [-]
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumePX(pressure float64, vapourFraction float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var v float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificVolumePX(p, vapourFraction)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of pressure for saturated liquid.
 *
 * @param pressure absolute pressure
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumeSaturatedLiquidP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var v float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificVolumeSaturatedLiquidP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of temperature for saturated liquid.
 *
 * @param temperature temperature
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumeSaturatedLiquidT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var v float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificVolumeSaturatedLiquidP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of pressure for saturated vapour.
 *
 * @param pressure absolute pressure
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumeSaturatedVapourP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var v float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificVolumeSaturatedVapourP(p)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of temperature for saturated vapour.
 *
 * @param temperature temperature
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumeSaturatedVapourT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var v float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificVolumeSaturatedVapourP(
		fourthRegion.REGION4.SaturationPressureT(T),
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Specific volume as a function of temperature &amp; vapour fraction.
 *
 * @param temperature temperature
 * @param vapourFraction vapour fraction [-]
 * @return specific volume
 * @return RangeError
 */
func (if97 *IF97) SpecificVolumeTX(temperature float64, vapourFraction float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var v float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	v = fourthRegion.REGION4.SpecificVolumePX(
		fourthRegion.REGION4.SaturationPressureT(T), vapourFraction,
	)

	return if97.ConvertFromDefault(if97.UnitSystem.SPECIFIC_VOLUME, v), err
}

/**
 * Speed of sound as a function of specific enthalpy &amp; specific entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return speed of sound
 * @return RangeError
 * @see #speedOfSoundPT(double, double)
 */
func (if97 *IF97) SpeedOfSoundHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var w float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}

	p := if97.region.PressureHS(h, s)
	T := if97.region.TemperatureHS(h, s)

	w = if97.region.SpeedOfSoundPT(p, T)
	return if97.ConvertFromDefault(if97.UnitSystem.SPEED_OF_SOUND, w), err
}

/**
 * Speed of sound as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return speed of sound
 * @return RangeError
 * @see #speedOfSoundPT(double, double)
 */
func (if97 *IF97) SpeedOfSoundPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var w float64

	reg, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}

	if reg == 4 {
		w = fourthRegion.REGION4.SpeedOfSoundPH(p, h)
	} else {
		T := if97.region.TemperaturePH(p, h)
		w = if97.region.SpeedOfSoundPT(p, T)
	}

	return if97.ConvertFromDefault(if97.UnitSystem.SPEED_OF_SOUND, w), err
}

/**
 * Speed of sound as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return speed of sound
 * @return RangeError
 * @see #speedOfSoundPT(double, double)
 */
func (if97 *IF97) SpeedOfSoundPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var w float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	T := if97.region.TemperaturePS(p, s)

	w = if97.region.SpeedOfSoundPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPEED_OF_SOUND, w), err
}

/**
 * Speed of sound as a function of pressure &amp; temperature.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return speed of sound
 * @return RangeError
 */
func (if97 *IF97) SpeedOfSoundPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var w float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	w = if97.region.SpeedOfSoundPT(p, T)

	return if97.ConvertFromDefault(if97.UnitSystem.SPEED_OF_SOUND, w), err
}

/**
 * Surface tension as a function of pressure.
 *
 * @param pressure absolute pressure
 * @return surface tension
 * @return RangeError
 */
func (if97 *IF97) SurfaceTensionP(pressure float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	var sigma float64

	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}

	sigma = fourthRegion.REGION4.SurfaceTensionT(
		fourthRegion.REGION4.SaturationTemperatureP(p),
	)
	return if97.ConvertFromDefault(if97.UnitSystem.SURFACE_TENSION, sigma), err
}

/**
 * Surface tension as a function of temperature.
 *
 * @param temperature temperature
 * @return surface tension
 * @return RangeError
 */
func (if97 *IF97) SurfaceTensionT(temperature float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var sigma float64

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}

	sigma = fourthRegion.REGION4.SurfaceTensionT(T)

	return if97.ConvertFromDefault(if97.UnitSystem.SURFACE_TENSION, sigma), err
}

/**
 * Temperature. [IF97 Supplementary Release S04]
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return temperature
 * @return RangeError
 */
func (if97 *IF97) TemperatureHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var T float64

	reg, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}
	T = if97.region.TemperatureHS(h, s)
	if reg==5 && math.IsNaN(T){
		return -1, errors.New("region 5 TemperatureHS unimplemented, Im sorry")
	}
	return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T), err
}

/**
 * Temperature.
 *
 * @param pressure absolute pressure [MPa]
 * @param enthalpy specific enthalpy [kJ/(kg)]
 * @return temperature [K]
 * @return RangeError
 */
func (if97 *IF97) TemperaturePH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var T float64

	reg, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}
	T = if97.region.TemperaturePH(p, h)
	if math.IsNaN(T) && reg==3{
		return -1, errors.New("unsupported subregion in region 3 TemperaturePH")
	}
	if math.IsNaN(T) && reg==5{
		return -1, errors.New("region 5 TemperaturePH unimplemented, Im sorry")
	}
	return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T), err
}

/**
 * Temperature.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return temperature
 * @return RangeError
 */
func (if97 *IF97) TemperaturePS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var T float64

	reg, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}
	T = if97.region.TemperaturePS(p, s)
	if math.IsNaN(T) && reg==3{
		return -1, errors.New("unsupported subregion in region 3 TemperaturePS")
	}
	if math.IsNaN(T) && reg==5{
		return -1, errors.New("TemperaturePS region 5 unimplemented, Im sorry")
	}

	
	return if97.ConvertFromDefault(if97.UnitSystem.TEMPERATURE, T), err
}

/**
 * Thermal conductivity as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return thermal conductivity
 * @return RangeError
 * @see #thermalConductivityRhoT(double, double)
 */
func (if97 *IF97) ThermalConductivityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var lambda float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}
	rho := 1.0 / if97.region.SpecificVolumeHS(h, s)
	T := if97.region.TemperatureHS(h, s)
	lambda = if97.defaultThermalConductivityRhoT(rho, T)

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_CONDUCTIVITY, lambda), err
}

/**
 * Thermal conductivity as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return thermal conductivity
 * @return RangeError
 * @see #thermalConductivityRhoT(double, double)
 */
func (if97 *IF97) ThermalConductivityPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var lambda float64

	_, err := if97.region.GetRegionPH(p, h)
	if err != nil {
		return -1, err
	}
	rho := 1.0 / if97.region.SpecificVolumePH(p, h)
	T := if97.region.TemperaturePH(p, h)
	lambda = if97.defaultThermalConductivityRhoT(rho, T)

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_CONDUCTIVITY, lambda), err
}

/**
 * Thermal conductivity as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return thermal conductivity
 * @return RangeError
 * @see #thermalConductivityRhoT(double, double)
 */
func (if97 *IF97) ThermalConductivityPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var lambda float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}
	rho := 1.0 / if97.region.SpecificVolumePS(p, s)
	T := if97.region.TemperaturePS(p, s)
	lambda = if97.defaultThermalConductivityRhoT(rho, T)

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_CONDUCTIVITY, lambda), err
}

/**
 * Thermal conductivity as a function of pressure &amp; temperature.
 *
 * Note that is method is not accurate in the two-phase region.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return thermal conductivity
 * @return RangeError
 * @see #thermalConductivityRhoT(double, double)
 */
func (if97 *IF97) ThermalConductivityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var lambda float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	rho := 1.0 / if97.region.SpecificVolumePT(p, T)
	lambda = if97.defaultThermalConductivityRhoT(rho, T)
	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_CONDUCTIVITY, lambda), err
}

/**
 * Thermal conductivity as a function of density &amp; temperature.
 *
 * @param density density
 * @param temperature temperature
 * @return thermal conductivity
 * @return RangeError
 */
func (if97 *IF97) ThermalConductivityRhoT(density float64, temperature float64) (float64, error) {

	rho := if97.ConvertToDefault(if97.UnitSystem.DENSITY, density)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var err error
	if T < constants.T0 {
		return -1, rangeError.ErrorFromValue(quantity.T, T, constants.T0)

	} else if T > constants.T5 {
		return -1, rangeError.ErrorFromValue(quantity.T, T, constants.T5)

	} else if rho <= 0 {
		return -1, rangeError.ErrorFromValue(quantity.Rho, rho, constants.T5)

	} else if rho > 1060 {
		return -1, rangeError.ErrorFromValue(quantity.Rho, rho, 1060)

	}

	var lambda float64

	lambda = if97.defaultThermalConductivityRhoT(rho, T)

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_CONDUCTIVITY, lambda), err
}

/**
 * Thermal diffusivity as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return thermal diffusivity
 * @return RangeError
 */
func (if97 *IF97) ThermalDiffusivityHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var kappa float64

	_, err := if97.region.GetRegionHS(h, s)
	if err != nil {
		return -1, err
	}
	p := if97.region.PressureHS(h, s)
	kappa, err = if97.defaultThermalDiffusivityPH(p, h)
	if err != nil {
		return -1, err
	}

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_DIFFUSIVITY, kappa), err
}

/**
 * Thermal diffusivity as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return thermal diffusivity
 * @return RangeError
 */
func (if97 *IF97) ThermalDiffusivityPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	var kappa float64

	kappa, err := if97.defaultThermalDiffusivityPH(p, h)
	if err != nil {
		return -1, err
	}

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_DIFFUSIVITY, kappa), err
}

/**
 * Thermal diffusivity as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return thermal diffusivity
 * @return RangeError
 */
func (if97 *IF97) ThermalDiffusivityPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)
	var kappa float64

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}
	h := if97.region.SpecificEnthalpyPS(p, s)
	kappa, err = if97.defaultThermalDiffusivityPH(p, h)
	if err != nil {
		return -1, err
	}

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_DIFFUSIVITY, kappa), err
}

/**
 * Thermal diffusivity as a function of pressure &amp; temperature.
 *
 * Note that is method is not accurate in the two-phase region.
 *
 * @param pressure absolute pressure
 * @param temperature temperature
 * @return thermal diffusivity
 * @return RangeError
 */
func (if97 *IF97) ThermalDiffusivityPT(pressure float64, temperature float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	var kappa float64

	_, err := if97.region.GetRegionPT(p, T)
	if err != nil {
		return -1, err
	}
	h := if97.region.SpecificEnthalpyPT(p, T)
	kappa, err = if97.defaultThermalDiffusivityPH(p, h)
	if err != nil {
		return -1, err
	}

	return if97.ConvertFromDefault(if97.UnitSystem.THERMAL_DIFFUSIVITY, kappa), err
}

/**
 * Vapour fraction as a function of specific enthalpy &amp; specific
 * entropy.
 *
 * @param enthalpy specific enthalpy
 * @param entropy specific entropy
 * @return vapour fraction [-]
 * @return RangeError
 */
func (if97 *IF97) VapourFractionHS(enthalpy float64, entropy float64) (float64, error) {

	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	err := fourthRegion.CheckHS(h, s)

	if err != nil {
		return -1, err
	}

	return fourthRegion.REGION4.VapourFractionHS(h, s), err

}

/**
 * Vapour fraction as a function of pressure &amp; specific enthalpy.
 *
 * @param pressure absolute pressure
 * @param enthalpy specific enthalpy
 * @return vapour fraction [-]
 * @return RangeError
 * @see #vapourFractionHS(double, double)
 */
func (if97 *IF97) VapourFractionPH(pressure float64, enthalpy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	h := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTHALPY, enthalpy)
	err := fourthRegion.CheckP(p)
	if err != nil {
		return -1, err
	}
	T := fourthRegion.REGION4.SaturationTemperatureP(p)
	err = fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}
	return fourthRegion.REGION4.VapourFractionPH(p, h), err
}

/**
 * Vapour fraction as a function of pressure &amp; specific entropy.
 *
 * @param pressure absolute pressure
 * @param entropy specific entropy
 * @return vapour fraction [-]
 * @return RangeError
 * @see #vapourFractionHS(double, double)
 */
func (if97 *IF97) VapourFractionPS(pressure float64, entropy float64) (float64, error) {

	p := if97.ConvertToDefault(if97.UnitSystem.PRESSURE, pressure)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	_, err := if97.region.GetRegionPS(p, s)
	if err != nil {
		return -1, err
	}

	return fourthRegion.REGION4.VapourFractionPS(p, s), err
}

/**
 * Vapour fraction as a function of temperature &amp; specific entropy.
 *
 * This method only returns values in the two-phase region.
 *
 * @param temperature temperature
 * @param entropy specific entropy
 * @return vapour fraction [-]
 * @return RangeError
 * @see #vapourFractionHS(double, double)
 */
func (if97 *IF97) VapourFractionTS(temperature float64, entropy float64) (float64, error) {

	T := if97.ConvertToDefault(if97.UnitSystem.TEMPERATURE, temperature)
	s := if97.ConvertToDefault(if97.UnitSystem.SPECIFIC_ENTROPY, entropy)

	err := fourthRegion.CheckT(T)
	if err != nil {
		return -1, err
	}
	ps := fourthRegion.REGION4.SaturationPressureT(T)

	err = fourthRegion.CheckP(ps)
	if err != nil {
		return -1, err
	}

	return fourthRegion.REGION4.VapourFractionTS(T, s), err

}
