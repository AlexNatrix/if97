package units

import (
	"fmt"
	"math"
	"strings"

	"if97.com/cmd/lib/utils/constants"
	"if97.com/cmd/lib/utils/quantity"
)

func setUnits(unitsList [22]string, units map[quantity.Quantity]string) map[quantity.Quantity]string {
	quantityList := [22]quantity.Quantity{
		quantity.P,
		quantity.T,
		quantity.V,
		quantity.Cp,
		quantity.Cv,
		quantity.Rho,
		quantity.Eta,
		quantity.Nu,
		quantity.A,
		quantity.U,
		quantity.H,
		quantity.S,
		quantity.Lambda,
		quantity.LambdaL,
		quantity.G,
		quantity.F,
		quantity.X,
		quantity.AlphaV,
		quantity.Kappa,
		quantity.KappaT,
		quantity.W,
		quantity.Sigma,
	}
	for i, v := range quantityList {
		units[v] = unitsList[i]
	}
	return units
}

//Default unit system.

var DEFAULT UnitSystem = UnitSystem{
	"DEFAULT",
	[]float64{1, 0}, // compressibility
	[]float64{1, 0}, // density
	[]float64{1, 0}, // dynamic viscosity
	[]float64{1, 0}, // isobaric cubic expansion coefficient
	[]float64{1, 0}, // kinematic viscosity
	[]float64{1, 0}, // pressure
	[]float64{1, 0}, // specific energy
	[]float64{1, 0}, // specific enthalpy
	[]float64{1, 0}, // specific entropy
	[]float64{1, 0}, // specific heat capacity
	[]float64{1, 0}, // specific volume
	[]float64{1, 0}, // speed of sound
	[]float64{1, 0}, // surface tension
	[]float64{1, 0}, // temperature
	[]float64{1, 0}, // thermal conductivity
	[]float64{1, 0}, // thermal diffusivity
	[]float64{1, 0}, // wavelength
}

var defaultUnits = [22]string{"MPa",
	"K",
	"m³/kg",
	"kJ/(kg·K)",
	"kJ/(kg·K)",
	"kg/m³",
	"Pa·s",
	"m²/s",
	"m²/s",
	"kJ/kg",
	"kJ/kg",
	"kJ/(kg·K)",
	"W/(m·K)",
	"μm",
	"kJ/kg",
	"kJ/kg",
	"",
	"1/K",
	"m²/s",
	"1/MPa",
	"m/s",
	"N/m",
}

// Engineering units.
var ENGINEERING UnitSystem = UnitSystem{
	"ENGINEERING",
	[]float64{1, 0},       // compressibility
	[]float64{1, 0},       // density
	[]float64{1, 0},       // dynamic viscosity
	[]float64{1, 0},       // isobaric cubic expansion coefficient
	[]float64{1, 0},       // kinematic viscosity
	[]float64{0.1, 0},     // pressure
	[]float64{1, 0},       // specific energy
	[]float64{1, 0},       // specific enthalpy
	[]float64{1, 0},       // specific entropy
	[]float64{1, 0},       // specific heat capacity
	[]float64{1, 0},       // specific volume
	[]float64{1, 0},       // speed of sound
	[]float64{1, 0},       // surface tension
	[]float64{1, constants.T0}, // temperature
	[]float64{1e3, 0},     // thermal conductivity
	[]float64{1, 0},       // thermal diffusivity
	[]float64{1, 0},       // wavelength
}

var engineeringUnits = [22]string{
	"bar",
	"°C",
	"m³/kg",
	"kJ/(kg·K)",
	"kJ/(kg·K)",
	"kg/m³",
	"Pa·s",
	"m²/s",
	"m²/s",
	"kJ/kg",
	"kJ/kg",
	"kJ/(kg·K)",
	"kW/(m·K)",
	"μm",
	"kJ/kg",
	"kJ/kg",
	"",
	"1/K",
	"m²/s",
	"1/MPa",
	"m/s",
	"N/m",
}

// SI unit system.
var SI UnitSystem = UnitSystem{
	"SI",
	[]float64{1e6, 0},  // compressibility
	[]float64{1, 0},    // density
	[]float64{1, 0},    // dynamic viscosity
	[]float64{1, 0},    // isobaric cubic expansion coefficient
	[]float64{1, 0},    // kinematic viscosity
	[]float64{1e-6, 0}, // pressure
	[]float64{1e-3, 0}, // specific energy
	[]float64{1e-3, 0}, // specific enthalpy
	[]float64{1e-3, 0}, // specific entropy
	[]float64{1e-3, 0}, // specific heat capacity
	[]float64{1, 0},    // specific volume
	[]float64{1, 0},    // speed of sound
	[]float64{1, 0},    // surface tension
	[]float64{1, 0},    // temperature
	[]float64{1, 0},    // thermal conductivity
	[]float64{1, 0},    // thermal diffusivity
	[]float64{1e6, 0},  // wavelength
}

var siUnits = [22]string{
	"Pa",
	"K",
	"m³/kg",
	"J/(kg·K)",
	"J/(kg·K)",
	"kg/m³",
	"Pa·s",
	"m²/s",
	"m²/s",
	"J/kg",
	"J/kg",
	"J/(kg·K)",
	"W/(m·K)",
	"m",
	"J/kg",
	"J/kg",
	"",
	"1/K",
	"m²/s",
	"1/Pa",
	"m/s",
	"N/m",
}

// British imperial unit system.
var IMPERIAL UnitSystem = UnitSystem{
	"IMPERIAL",
	[]float64{1e6 * constants.InSquare / constants.Lb, 0}, // compressibility
	[]float64{constants.Lb / constants.FtCube, 0},         // density
	[]float64{1e-3, 0},                                                               // dynamic viscosity
	[]float64{1 / constants.Ra, 0},                                                   // isobaric cubic expansion coefficient
	[]float64{1e-6, 0},                                                               // kinematic viscosity [centiStokes]
	[]float64{constants.Psi, 0},                                                      // pressure
	[]float64{constants.BTU / constants.Lb, 0},                                       // specific energy
	[]float64{constants.BTU / constants.Lb, 0},                                       // specific enthalpy
	[]float64{constants.BTU / (constants.Lb * constants.Ra), 0},                      // specific entropy
	[]float64{constants.BTU / (constants.Lb * constants.Ra), 0},                      // specific heat capacity
	[]float64{constants.FtCube / constants.Lb, 0},                                    // specific volume
	[]float64{constants.Ft, 0},                                                       // speed of sound
	[]float64{constants.Lbf / constants.Ft, 0},                                       // surface tension
	[]float64{5 / 9, 459.67 * 5 / 9},                                                 // temperature
	[]float64{1e3 * constants.BTU / (constants.Hr * constants.Ft * constants.Ra), 0}, // thermal conductivity
	[]float64{1e-6, 0},                                                               // thermal diffusivity
	[]float64{constants.In * 1e6, 0},                                                 // wavelength
}

var imperialUnits = [22]string{
	"psi",
	"°F",
	"ft³/lb",
	"BTU/(lb·R)",
	"BTU/(lb·R)",
	"lb/ft³",
	"cP",
	"cSt",
	"cSt",
	"BTU/lb",
	"BTU/lb",
	"BTU/(lb·R)",
	"BTU/(hr·ft·R)",
	"in",
	"BTU/lb",
	"BTU/lb",
	"",
	"1/R",
	"cSt",
	"in²/lb",
	"ft/s",
	"lbf/ft",
}

type UnitSystem struct {
	name                                 string
	COMPRESSIBILITY                      []float64
	DENSITY                              []float64
	DYNAMIC_VISCOSITY                    []float64
	ISOBARIC_CUBIC_EXPANSION_COEFFICIENT []float64
	KINEMATIC_VISCOSITY                  []float64
	PRESSURE                             []float64
	SPECIFIC_ENERGY                      []float64
	SPECIFIC_ENTHALPY                    []float64
	SPECIFIC_ENTROPY                     []float64
	SPECIFIC_HEAT_CAPACITY               []float64
	SPECIFIC_VOLUME                      []float64
	SPEED_OF_SOUND                       []float64
	SURFACE_TENSION                      []float64
	TEMPERATURE                          []float64
	THERMAL_CONDUCTIVITY                 []float64
	THERMAL_DIFFUSIVITY                  []float64
	WAVELENGTH                           []float64
}

var UNITS map[quantity.Quantity]string

var tol = 1e-9

func (us *UnitSystem) apply() {
	if math.Abs(us.PRESSURE[0]-1) < tol {
		/*
		 Default
		*/
		setUnits(defaultUnits, UNITS)

	} else if math.Abs(us.PRESSURE[0]-0.1) < tol {
		/*
		 Engineering
		*/
		setUnits(engineeringUnits, UNITS)

	} else if math.Abs(us.PRESSURE[0]-1e-6) < tol {
		/*
		 SI
		*/
		setUnits(siUnits, UNITS)
	} else if math.Abs(us.PRESSURE[0]-constants.Psi) < tol {
		/*
		 Imperial
		*/
		setUnits(imperialUnits, UNITS)
	} else {
		fmt.Errorf("Unsupported unit system")
		
	}
}

func getLabel(quantity quantity.Quantity) string {

	if _, ok := UNITS[quantity]; !ok {
		fmt.Errorf("Unknown label for quantity: " + quantity.String())
		return fmt.Sprintf("Unknown label for quantity %s", quantity.String())
	}
	return fmt.Sprintf("%s [%s]", quantity.String(), getUnit(quantity))
}

func getUnit(quantity quantity.Quantity) string {
	if _, ok := UNITS[quantity]; ok {
		return UNITS[quantity]
	}
	return "-"
}

func (us *UnitSystem) String() string {
	switch us {
	case &SI:
		return us.name

	default:
		return strings.ToLower(us.name)
	}
}
