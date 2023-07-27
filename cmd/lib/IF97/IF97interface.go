package IF97

import (
	"if97.com/cmd/lib/utils/quantity"
	"if97.com/cmd/lib/utils/units"
)

/**
 * <p>
 * Steam tables for industrial use according to the international standard for
 * the properties of water and steam, the IAPWS-IF97 formulation and the
 * international standards for transport and other properties.</p>
 *
 * <p>
 * By default, units are as given in the book cited below. See
 * {@link UnitSystem} for options.</p>
 *
 * <p>
 * Properties are generally available as functions of three quantity
 * combinations: pressure &amp; temperature (p, T), pressure &amp; specific
 * enthalpy (p, h), and specific enthalpy &amp; specific entropy (h, s).</p>
 *
 * <ul> <li>Wagner, Wolfgang &amp; Kretzschmar, Hans-Joachim, 2008,
 * <i>International Steam Tables &mdash; Properties of Water and Steam Based on
 * the Industrial Formulation IAPWS-IF97</i>, 2<sup>nd</sup> Edition,
 * Springer-Verlag, Berlin Heidelberg, ISBN 978-3-540-21419-9</li> </ul>
 *
 */
type IF97calc interface {
	/**
	 * Prandtl number.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return Prandtl number
	 * @throws OutOfRangeException out-of-range exception
	 */
	PrandtlHS(enthalpy float64, entropy float64) (float64, error)

	/**
	 * Prandtl number.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return Prandtl number
	 * @throws OutOfRangeException out-of-range exception
	 */
	PrandtlPH(pressure float64, enthalpy float64) (float64, error)

	/**
	 * Prandtl number.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return Prandtl number
	 * @throws OutOfRangeException out-of-range exception
	 */
	PrandtlPS(pressure float64, entropy float64) (float64, error)

	/**
	 * Prandtl number.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return Prandtl number
	 * @throws OutOfRangeException out-of-range exception
	 */
	PrandtlPT(pressure float64, temperature float64) (float64, error)

	/**
	 * Isothermal compressibility as a function of specific enthalpy &amp;
	 * specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return isothermal compressibility
	 * @throws OutOfRangeException out-of-range exception
	 */
	CompressibilityHS(enthalpy float64, entropy float64) (float64, error)

	/**
	 * Isothermal compressibility as a function of pressure &amp; specific
	 * enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return isothermal compressibility
	 * @throws OutOfRangeException out-of-range exception
	 */
	CompressibilityPH(pressure float64, enthalpy float64) (float64, error)

	/**
	 * Isothermal compressibility as a function of pressure &amp; specific
	 * entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return isothermal compressibility
	 * @throws OutOfRangeException out-of-range exception
	 */
	CompressibilityPS(pressure float64, entropy float64) (float64, error)

	/**
	 * Isothermal compressibility as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return isothermal compressibility
	 * @throws OutOfRangeException out-of-range exception
	 */
	CompressibilityPT(pressure float64, temperature float64) (float64, error)

	/**
	 * Compression factor (real-gas factor) as a function of pressure &amp;
	 * temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return compression factor
	 * @throws OutOfRangeException out-of-range exception
	 */
	// CompressionFactorPT(pressure float64, temperature float64) (float64, error)

	// ConvertFromDefaultQuantity(us units.UnitSystem, quantity quantity.Quantity, value float64)

	// ConvertFromDefault(quantity []float64, value float64)

	// ConvertToDefault(quantity []float64, value float64)

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
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificVolumeHS(double, double)
	 */
	DensityHS(enthalpy float64, entropy float64) (float64, error)
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
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificVolumePH(double, double)
	 */
	DensityPH(pressure float64, enthalpy float64) (float64, error)

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
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificVolumePS(double, double)
	 */
	DensityPS(pressure float64, entropy float64) (float64, error)
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
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificVolumePT(double, double)
	 */
	DensityPT(pressure float64, temperature float64) (float64, error)
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
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificVolumePX(double, double)
	 */
	DensityPX(pressure float64, vapourFraction float64) (float64, error)
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
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificVolumeTX(double, double)
	 */
	DensityTX(temperature float64, vapourFraction float64) (float64, error)
	/**
	 * Dielectric constant as a function of specific enthalpy and specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return dielectric constant [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dielectricConstantRhoT(double, double)
	 */
	DielectricConstantHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Dielectric constant as a function of pressure and specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return dielectric constant [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dielectricConstantRhoT(double, double)
	 */
	DielectricConstantPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Dielectric constant as a function of pressure and specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return dielectric constant [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dielectricConstantRhoT(double, double)
	 */
	DielectricConstantPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Dielectric constant as a function of pressure and temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return dielectric constant [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dielectricConstantRhoT(double, double)
	 */
	DielectricConstantPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Dielectric constant (relative static dielectric constant or relative
	 * static permittivity) as a function of density and temperature.
	 *
	 * @param density density
	 * @param temperature temperature
	 * @return dielectric constant [-]
	 * @throws OutOfRangeException out-of-range exception
	 */
	DielectricConstantRhoT(density float64, temperature float64) (float64, error)
	/**
	 * Dynamic viscosity as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return dynamic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dynamicViscosityRhoT(double, double)
	 */
	DynamicViscosityHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Dynamic viscosity as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return dynamic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dynamicViscosityRhoT(double, double)
	 */
	DynamicViscosityPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Dynamic viscosity as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return dynamic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dynamicViscosityRhoT(double, double)
	 */
	DynamicViscosityPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Dynamic viscosity as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return dynamic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #dynamicViscosityRhoT(double, double)
	 */
	DynamicViscosityPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Dynamic viscosity as a function of density &amp; temperature.
	 *
	 * @param density density
	 * @param temperature temperature
	 * @return dynamic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 */
	DynamicViscosityRhoT(density float64, temperature float64) (float64, error)

	getRegionPT(pressure float64, temperature float64) (int, error)
	/**
	 * Gets the unit system.
	 *
	 * @return unit system
	 */
	getUnitSystem() units.UnitSystem
	/**
	 * Heat capacity ratio as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return heat capacity ratio
	 * @throws OutOfRangeException out-of-range exception
	 * @see #heatCapacityRatioPT(double, double)
	 */
	HeatCapacityRatioHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Heat capacity ratio as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return heat capacity ratio
	 * @throws OutOfRangeException out-of-range exception
	 * @see #heatCapacityRatioPT(double, double)
	 */
	HeatCapacityRatioPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Heat capacity ratio as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return heat capacity ratio
	 * @throws OutOfRangeException out-of-range exception
	 * @see #heatCapacityRatioPT(double, double)
	 */
	HeatCapacityRatioPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Heat capacity ratio as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return heat capacity ratio
	 * @throws OutOfRangeException out-of-range exception
	 */
	HeatCapacityRatioPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Isentropic exponent as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return isentropic exponent
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isentropicExponentPT(double, double)
	 */
	IsentropicExponentHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Isentropic exponent as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return isentropic exponent
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isentropicExponentPT(double, double)
	 */
	IsentropicExponentPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Isentropic exponent as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return isentropic exponent
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isentropicExponentPT(double, double)
	 */
	IsentropicExponentPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Isentropic exponent as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return isentropic exponent
	 * @throws OutOfRangeException out-of-range exception
	 */
	IsentropicExponentPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Isobaric cubic expansion coefficient as a function of specific enthalpy
	 * &amp; specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return isobaric cubic expansion coefficient
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricCubicExpansionCoefficientPT(double, double)
	 */
	IsobaricCubicExpansionCoefficientHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Isobaric cubic expansion coefficient as a function of pressure &amp;
	 * specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return isobaric cubic expansion coefficient
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricCubicExpansionCoefficientPT(double, double)
	 */
	IsobaricCubicExpansionCoefficientPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Isobaric cubic expansion coefficient as a function of pressure &amp;
	 * specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return isobaric cubic expansion coefficient
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricCubicExpansionCoefficientPT(double, double)
	 */
	IsobaricCubicExpansionCoefficientPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Isobaric cubic expansion coefficient as a function of pressure &amp;
	 * temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return isobaric cubic expansion coefficient
	 * @throws OutOfRangeException out-of-range exception
	 */
	IsobaricCubicExpansionCoefficientPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Isobaric cubic expansion coefficient as a function of pressure &amp;
	 * vapour fraction.
	 *
	 * @param pressure absolute pressure
	 * @param vapourFraction vapour fraction
	 * @return isobaric cubic expansion coefficient
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricCubicExpansionCoefficientPH(double, double)
	 */
	IsobaricCubicExpansionCoefficientPX(pressure float64, vapourFraction float64)
	/**
	 * Isobaric cubic expansion coefficient as a function of temperature &amp;
	 * vapour fraction.
	 *
	 * @param temperature temperature
	 * @param vapourFraction vapour fraction
	 * @return isobaric cubic expansion coefficient
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricCubicExpansionCoefficientPH(double, double)
	 */
	IsobaricCubicExpansionCoefficientTX(temperature float64, vapourFraction float64)
	/**
	 * Specific isobaric heat capacity as a function of specific enthalpy &amp;
	 * specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return specific isobaric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricHeatCapacityPT(double, double)
	 */
	IsobaricHeatCapacityHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Specific isobaric heat capacity as a function of pressure &amp; specific
	 * enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return specific isobaric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricHeatCapacityPT(double, double)
	 */
	IsobaricHeatCapacityPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Specific isobaric heat capacity as a function of pressure &amp; specific
	 * entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return specific isobaric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isobaricHeatCapacityPT(double, double)
	 */
	IsobaricHeatCapacityPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Specific isobaric heat capacity as a function of pressure &amp;
	 * temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific isobaric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 */
	IsobaricHeatCapacityPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Specific isochoric heat capacity as a function of specific enthalpy &amp;
	 * specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return specific isochoric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isochoricHeatCapacityPT(double, double)
	 */
	IsochoricHeatCapacityHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Specific isochoric heat capacity as a function of pressure &amp; specific
	 * enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return specific isochoric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isochoricHeatCapacityPT(double, double)
	 */
	IsochoricHeatCapacityPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Specific isochoric heat capacity as a function of pressure &amp; specific
	 * entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return specific isochoric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #isochoricHeatCapacityPT(double, double)
	 */
	IsochoricHeatCapacityPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Specific isochoric heat capacity as a function of pressure &amp;
	 * temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific isochoric heat capacity
	 * @throws OutOfRangeException out-of-range exception
	 */
	IsochoricHeatCapacityPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Kinematic viscosity as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return kinematic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #kinematicViscosityRhoT(double, double)
	 */
	KinematicViscosityHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Kinematic viscosity as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return kinematic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #kinematicViscosityRhoT(double, double)
	 */
	KinematicViscosityPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Kinematic viscosity as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return kinematic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #kinematicViscosityRhoT(double, double)
	 */
	KinematicViscosityPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Kinematic viscosity as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return kinematic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #kinematicViscosityRhoT(double, double)
	 */
	KinematicViscosityPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Kinematic viscosity as a function of density &amp; temperature.
	 *
	 * @param density density
	 * @param temperature temperature
	 * @return kinematic viscosity
	 * @throws OutOfRangeException out-of-range exception
	 */
	KinematicViscosityRhoT(density float64, temperature float64) (float64, error)
	/**
	 * Partial derivative of z with respect to x for constant y in SI units, as
	 * a function of pressure and specific enthalpy.
	 *
	 * <p>
	 * (<sup>&part;z</sup>/<sub>&part;x</sub>)<sub>y</sub>(p, h)
	 * </p>
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @param x any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @param y any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @param z any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @return partial derivative [SI units]
	 * @throws OutOfRangeException out-of-range exception
	 * @see Quantity#getPartialDerivatives()
	 */
	PartialDerivativePH(pressure float64, enthalpy float64, x quantity.Quantity, y quantity.Quantity, z quantity.Quantity) (float64, error)
	/**
	 * Gets the partial derivative of z with respect to x for constant y in SI
	 * units, as a function of pressure and temperature.
	 *
	 * <p>
	 * (<sup>&part;z</sup>/<sub>&part;x</sub>)<sub>y</sub>(p, T)
	 * </p>
	 *
	 * Note that this method is not suitable for the saturated region as
	 * pressure and temperature are coupled. Preferably, use
	 * {@link #partialDerivativePH(double, double, com.hummeling.if97.IF97.Quantity, com.hummeling.if97.IF97.Quantity, com.hummeling.if97.IF97.Quantity)}
	 * instead.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @param x any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @param y any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @param z any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @return partial derivative [SI units]
	 * @throws OutOfRangeException out-of-range exception
	 * @see Quantity#getPartialDerivatives()
	 * @see #partialDerivativePH(double, double,
	 */
	PartialDerivativePT(pressure float64, temperature float64, x quantity.Quantity, y quantity.Quantity, z quantity.Quantity) (float64, error)
	/**
	 * Gets the partial derivative of z with respect to x for constant y in SI
	 * units, as a function of density and temperature, valid in region 3 only!
	 *
	 * <p>
	 * (<sup>&part;z</sup>/<sub>&part;x</sub>)<sub>y</sub>(&rho;, T)
	 * </p>
	 * Preferably, use
	 * {@link #partialDerivativePH(double, double, com.hummeling.if97.IF97.Quantity, com.hummeling.if97.IF97.Quantity, com.hummeling.if97.IF97.Quantity)}
	 * instead.
	 *
	 * @param density density
	 * @param temperature temperature
	 * @param x any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @param y any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @param z any {@link Quantity} part of the set returned by
	 * {@link Quantity#getPartialDerivatives()}
	 * @return partial derivative [SI units]
	 * @throws OutOfRangeException out-of-range exception
	 * @see Quantity#getPartialDerivatives()
	 * @see #partialDerivativePH(double, double,
	 * com.hummeling.if97.IF97.Quantity, com.hummeling.if97.IF97.Quantity,
	 * com.hummeling.if97.IF97.Quantity)
	 */
	PartialDerivativeRhoT(density float64, temperature float64, x quantity.Quantity, y quantity.Quantity, z quantity.Quantity) (float64, error)

	/**
	 * Pressure as a function of specific enthalpy &amp; specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return pressure
	 * @throws OutOfRangeException out-of-range exception
	 */
	PressureHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Refractive index as a function of specific enthalpy, specific entropy
	 * &amp; wave length.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @param wavelength wavelength
	 * @return refractive index [-]
	 * @throws OutOfRangeException out-of-range exception
	 */
	RefractiveIndexHSLambda(enthalpy float64, entropy float64, wavelength float64) (float64, error)
	/**
	 * Refractive index as a function of pressure, specific enthalpy &amp; wave
	 * length.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @param wavelength wavelength
	 * @return refractive index [-]
	 * @throws OutOfRangeException out-of-range exception
	 */
	RefractiveIndexPHLambda(pressure float64, enthalpy float64, wavelength float64) (float64, error)
	/**
	 * Refractive index as a function of pressure, specific entropy &amp; wave
	 * length.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @param wavelength wavelength
	 * @return refractive index [-]
	 * @throws OutOfRangeException out-of-range exception
	 */
	RefractiveIndexPSLambda(pressure float64, entropy float64, wavelength float64) (float64, error)
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
	RefractiveIndexPTLambda(pressure float64, temperature float64, wavelength float64) (float64, error)

	/**
	 * Refractive index as a function of density, temperature &amp; wave length.
	 *
	 * @param density density
	 * @param temperature temperature
	 * @param waveLength wave length
	 * @return refractive index [-]
	 * @throws OutOfRangeException out-of-range exception
	 */
	RefractiveIndexRhoTLambda(density float64, temperature float64, waveLength float64) (float64, error)
	/**
	 * Saturation pressure as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return saturation pressure
	 * @throws OutOfRangeException out-of-range exception
	 */
	SaturationPressureHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Saturation pressure as a function of temperature.
	 *
	 * @param temperature saturation temperature
	 * @return saturation pressure
	 * @throws OutOfRangeException out-of-range exception
	 */
	SaturationPressureT(temperature float64) (float64, error)
	/**
	 * Saturation temperature as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return saturation temperature
	 * @throws OutOfRangeException out-of-range exception
	 */
	SaturationTemperatureHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Saturation temperature as a function of pressure.
	 *
	 * @param pressure saturation pressure
	 * @return saturation temperature
	 * @throws OutOfRangeException out-of-range exception
	 */
	SaturationTemperatureP(pressure float64) (float64, error)
	/**
	 * Sets (changes) the unit system.
	 *
	 * @param unitSystem unit system
	 */
	setUnitSystem(unitSystem units.UnitSystem)
	/**
	 * Specific enthalpy as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEnthalpyPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Specific enthalpy as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEnthalpyPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Specific enthalpy as a function of pressure &amp; vapour fraction.
	 *
	 * @param pressure absolute pressure
	 * @param vapourFraction vapour fraction [-]
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEnthalpyPX(pressure float64, vapourFraction float64) (float64, error)
	/**
	 * Specific enthalpy as a function of pressure for saturated liquid.
	 *
	 * @param pressure saturation pressure
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEnthalpyPX(double, double)
	 */
	SpecificEnthalpySaturatedLiquidP(pressure float64) (float64, error)
	/**
	 * Specific enthalpy as a function of temperature for saturated liquid.
	 *
	 * @param temperature saturation temperature
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEnthalpyTX(double, double)
	 */
	SpecificEnthalpySaturatedLiquidT(temperature float64) (float64, error)
	/**
	 * Specific enthalpy as a function of pressure for saturated vapour.
	 *
	 * @param pressure saturation pressure
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEnthalpyPX(double, double)
	 */
	SpecificEnthalpySaturatedVapourP(pressure float64) (float64, error)
	/**
	 * Specific enthalpy as a function of temperature for saturated vapour.
	 *
	 * @param temperature saturation temperature
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEnthalpyTX(double, double)
	 */
	SpecificEnthalpySaturatedVapourT(temperature float64) (float64, error)
	/**
	 * Specific enthalpy as a function of temperature &amp; vapour fraction.
	 *
	 * @param temperature temperature
	 * @param vapourFraction vapour fraction [-]
	 * @return specific enthalpy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEnthalpyTX(temperature float64, vapourFraction float64) (float64, error)
	/**
	 * Specific entropy as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEntropyPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Specific entropy as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEntropyPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Specific entropy as a function of pressure &amp; vapour fraction.
	 *
	 * @param pressure absolute pressure
	 * @param vapourFraction vapour fraction [-]
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEntropyPX(pressure float64, vapourFraction float64) (float64, error)
	/**
	 * Specific entropy as a function of pressure for saturated liquid.
	 *
	 * @param pressure saturation pressure
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEntropyPX(double, double)
	 */
	SpecificEntropySaturatedLiquidP(pressure float64) (float64, error)
	/**
	 * Specific entropy as a function of temperature for saturated liquid.
	 *
	 * @param temperature saturation temperature
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEntropyTX(double, double)
	 */
	SpecificEntropySaturatedLiquidT(temperature float64) (float64, error)
	/**
	 * Specific entropy as a function of pressure for saturated vapour.
	 *
	 * @param pressure saturation pressure
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEntropyPX(double, double)
	 */
	SpecificEntropySaturatedVapour(pressure float64) (float64, error)
	/**
	 * Specific entropy as a function of temperature for saturated vapour.
	 *
	 * @param temperature saturation temperature
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificEntropyTX(double, double)
	 */
	SpecificEntropySaturatedVapourT(temperature float64) (float64, error)
	/**
	 * Specific entropy as a function of temperature &amp; vapour fraction.
	 *
	 * @param temperature temperature
	 * @param vapourFraction vapour fraction [-]
	 * @return specific entropy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificEntropyTX(temperature float64, vapourFraction float64) (float64, error)
	/**
	 * Specific Gibbs free energy as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific Gibbs free energy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificGibbsFreeEnergyPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Specific internal energy as a function of specific enthalpy &amp;
	 * specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyPT(double, double)
	 */
	SpecificInternalEnergyHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Specific internal energy as a function of pressure &amp; specific
	 * enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyPT(double, double)
	 */
	SpecificInternalEnergyPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Specific internal energy as a function of pressure &amp; specific
	 * entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyPT(double, double)
	 */
	SpecificInternalEnergyPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Specific internal energy as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificInternalEnergyPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Specific internal energy as a function of pressure &amp; vapour fraction.
	 *
	 * @param pressure absolute pressure
	 * @param vapourFraction vapour fraction [-]
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificInternalEnergyPX(pressure float64, vapourFraction float64) (float64, error)
	/**
	 * Specific internal energy as a function of pressure for saturated liquid.
	 *
	 * @param pressure saturation pressure
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyPX(double, double)
	 */
	SpecificInternalEnergySaturatedLiquidP(pressure float64) (float64, error)
	/**
	 * Specific internal energy as a function of temperature for saturated
	 * liquid.
	 *
	 * @param temperature saturation temperature
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyTX(double, double)
	 */
	SpecificInternalEnergySaturatedLiquidT(temperature float64) (float64, error)
	/**
	 * Specific internal energy as a function of pressure for saturated vapour.
	 *
	 * @param pressure saturation pressure
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyPX(double, double)
	 */
	SpecificInternalEnergySaturatedVapourP(pressure float64) (float64, error)
	/**
	 * Specific internal energy as a function of temperature for saturated
	 * vapour.
	 *
	 * @param temperature saturation temperature
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 * @see #specificInternalEnergyTX(double, double)
	 */
	SpecificInternalEnergySaturatedVapourT(temperature float64) (float64, error)
	/**
	 * Specific internal energy as a function of temperature &amp; vapour
	 * fraction.
	 *
	 * @param temperature temperature
	 * @param vapourFraction vapour fraction [-]
	 * @return specific internal energy
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificInternalEnergyTX(temperature float64, vapourFraction float64) (float64, error)
	/**
	 * Specific volume as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumeHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Specific volume as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumePH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Specific volume as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumePS(pressure float64, entropy float64) (float64, error)
	/**
	 * Specific volume as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumePT(pressure float64, temperature float64) (float64, error)
	/**
	 * Specific volume as a function of pressure &amp; vapour fraction.
	 *
	 * @param pressure absolute pressure
	 * @param vapourFraction vapour fraction [-]
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumePX(pressure float64, vapourFraction float64) (float64, error)
	/**
	 * Specific volume as a function of pressure for saturated liquid.
	 *
	 * @param pressure absolute pressure
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumeSaturatedLiquidP(pressure float64) (float64, error)
	/**
	 * Specific volume as a function of temperature for saturated liquid.
	 *
	 * @param temperature temperature
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumeSaturatedLiquidT(temperature float64) (float64, error)
	/**
	 * Specific volume as a function of pressure for saturated vapour.
	 *
	 * @param pressure absolute pressure
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumeSaturatedVapourP(pressure float64) (float64, error)
	/**
	 * Specific volume as a function of temperature for saturated vapour.
	 *
	 * @param temperature temperature
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumeSaturatedVapourT(temperature float64) (float64, error)
	/**
	 * Specific volume as a function of temperature &amp; vapour fraction.
	 *
	 * @param temperature temperature
	 * @param vapourFraction vapour fraction [-]
	 * @return specific volume
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpecificVolumeTX(temperature float64, vapourFraction float64) (float64, error)
	/**
	 * Speed of sound as a function of specific enthalpy &amp; specific entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return speed of sound
	 * @throws OutOfRangeException out-of-range exception
	 * @see #speedOfSoundPT(double, double)
	 */
	SpeedOfSoundHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Speed of sound as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return speed of sound
	 * @throws OutOfRangeException out-of-range exception
	 * @see #speedOfSoundPT(double, double)
	 */
	SpeedOfSoundPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Speed of sound as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return speed of sound
	 * @throws OutOfRangeException out-of-range exception
	 * @see #speedOfSoundPT(double, double)
	 */
	SpeedOfSoundPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Speed of sound as a function of pressure &amp; temperature.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return speed of sound
	 * @throws OutOfRangeException out-of-range exception
	 */
	SpeedOfSoundPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Surface tension as a function of pressure.
	 *
	 * @param pressure absolute pressure
	 * @return surface tension
	 * @throws OutOfRangeException out-of-range exception
	 */
	SurfaceTensionP(pressure float64) (float64, error)
	/**
	 * Surface tension as a function of temperature.
	 *
	 * @param temperature temperature
	 * @return surface tension
	 * @throws OutOfRangeException out-of-range exception
	 */
	SurfaceTensionT(temperature float64) (float64, error)
	/**
	 * Temperature. [IF97 Supplementary Release S04]
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return temperature
	 * @throws OutOfRangeException out-of-range exception
	 */
	TemperatureHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Temperature.
	 *
	 * @param pressure absolute pressure [MPa]
	 * @param enthalpy specific enthalpy [kJ/(kg)]
	 * @return temperature [K]
	 * @throws OutOfRangeException out-of-range exception
	 */
	TemperaturePH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Temperature.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return temperature
	 * @throws OutOfRangeException out-of-range exception
	 */
	TemperaturePS(pressure float64, entropy float64) (float64, error)
	/**
	 * Thermal conductivity as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return thermal conductivity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #thermalConductivityRhoT(double, double)
	 */
	ThermalConductivityHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Thermal conductivity as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return thermal conductivity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #thermalConductivityRhoT(double, double)
	 */
	ThermalConductivityPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Thermal conductivity as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return thermal conductivity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #thermalConductivityRhoT(double, double)
	 */
	ThermalConductivityPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Thermal conductivity as a function of pressure &amp; temperature.
	 *
	 * Note that is method is not accurate in the two-phase region.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return thermal conductivity
	 * @throws OutOfRangeException out-of-range exception
	 * @see #thermalConductivityRhoT(double, double)
	 */
	ThermalConductivityPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Thermal conductivity as a function of density &amp; temperature.
	 *
	 * @param density density
	 * @param temperature temperature
	 * @return thermal conductivity
	 * @throws OutOfRangeException out-of-range exception
	 */
	ThermalConductivityRhoT(density float64, temperature float64) (float64, error)
	/**
	 * Thermal diffusivity as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return thermal diffusivity
	 * @throws OutOfRangeException out-of-range exception
	 */
	ThermalDiffusivityHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Thermal diffusivity as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return thermal diffusivity
	 * @throws OutOfRangeException out-of-range exception
	 */
	ThermalDiffusivityPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Thermal diffusivity as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return thermal diffusivity
	 * @throws OutOfRangeException out-of-range exception
	 */
	ThermalDiffusivityPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Thermal diffusivity as a function of pressure &amp; temperature.
	 *
	 * Note that is method is not accurate in the two-phase region.
	 *
	 * @param pressure absolute pressure
	 * @param temperature temperature
	 * @return thermal diffusivity
	 * @throws OutOfRangeException out-of-range exception
	 */
	ThermalDiffusivityPT(pressure float64, temperature float64) (float64, error)
	/**
	 * Vapour fraction as a function of specific enthalpy &amp; specific
	 * entropy.
	 *
	 * @param enthalpy specific enthalpy
	 * @param entropy specific entropy
	 * @return vapour fraction [-]
	 * @throws OutOfRangeException out-of-range exception
	 */
	VapourFractionHS(enthalpy float64, entropy float64) (float64, error)
	/**
	 * Vapour fraction as a function of pressure &amp; specific enthalpy.
	 *
	 * @param pressure absolute pressure
	 * @param enthalpy specific enthalpy
	 * @return vapour fraction [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #vapourFractionHS(double, double)
	 */
	VapourFractionPH(pressure float64, enthalpy float64) (float64, error)
	/**
	 * Vapour fraction as a function of pressure &amp; specific entropy.
	 *
	 * @param pressure absolute pressure
	 * @param entropy specific entropy
	 * @return vapour fraction [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #vapourFractionHS(double, double)
	 */
	VapourFractionPS(pressure float64, entropy float64) (float64, error)
	/**
	 * Vapour fraction as a function of temperature &amp; specific entropy.
	 *
	 * This method only returns values in the two-phase region.
	 *
	 * @param temperature temperature
	 * @param entropy specific entropy
	 * @return vapour fraction [-]
	 * @throws OutOfRangeException out-of-range exception
	 * @see #vapourFractionHS(double, double)
	 */
	VapourFractionTS(temperature float64, entropy float64) (float64, error)
}
