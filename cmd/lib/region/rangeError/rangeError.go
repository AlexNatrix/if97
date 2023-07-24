package rangeError

import (
	"errors"
	"fmt"
	"strings"

	"if97.com/cmd/lib/utils/quantity"
	"if97.com/cmd/lib/utils/units"
)

type RangeError struct {
	QUANTITIES  []quantity.Quantity
	VALUES      []float64
	LIMITS      []float64
	UNIT_SYSTEM units.UnitSystem
}


// 	/**
// 	 * Get exceeded limit.
// 	 *
// 	 * @return limit value
// 	 */
func (err *RangeError) GetLimit() float64 {
	return err.LIMITS[0]
}


// 	/**
// 	 * Get exceeded Quantity.
// 	 *
// 	 * @return Quantity value
// 	 */
func (err *RangeError) GetQuantity() string {
		return err.QUANTITIES[0].String();
}


func ConvertFromDefault(quantity []float64, value float64)float64 {
	return (value - quantity[1]) / quantity[0];
}

func ConvertFromDefaultQuantity(unitSystem units.UnitSystem, quantity quantity.Quantity, value float64) (float64, error) {

	switch quantity.NAME {
	case "temperature":
		return ConvertFromDefault(unitSystem.TEMPERATURE, value),nil

	case "specific Helmholtz free energy":
	case "specific Gibbs free energy":
	case "specific internal energy":
		return ConvertFromDefault(unitSystem.SPECIFIC_ENERGY, value),nil

	case "specific enthalpy":
		return ConvertFromDefault(unitSystem.SPECIFIC_ENTHALPY, value),nil

	case "thermal conductivity":
		return ConvertFromDefault(unitSystem.THERMAL_CONDUCTIVITY, value),nil

	case "wavelength":
		return ConvertFromDefault(unitSystem.WAVELENGTH, value),nil

	case "absolute pressure":
		return ConvertFromDefault(unitSystem.PRESSURE, value),nil

	case "density":
		return ConvertFromDefault(unitSystem.DENSITY, value),nil

	case "specific entropy":
		return ConvertFromDefault(unitSystem.SPECIFIC_ENTROPY, value),nil

	case "specific volume":
		return ConvertFromDefault(unitSystem.SPECIFIC_VOLUME, value),nil

	case "vapour fraction":
		return value,nil

	}
	return -1, errors.New("no conversion available for: " + quantity.String())
}

func  ErrorFromValue(q quantity.Quantity , value float64, limit float64) RangeError{
	return RangeError{
		[]quantity.Quantity{q}, 
		[]float64{value}, 
		[]float64{limit},
		units.DEFAULT,
		};
}


func (re *RangeError) ConvertFromDefault(unitSystem units.UnitSystem ) error {

	var values []float64
	var limits  []float64
	var err error
for i := 0; i < len(re.QUANTITIES); i++ {
	values[i],err = ConvertFromDefaultQuantity(unitSystem  , re.QUANTITIES[i], re.VALUES[i]);
	if err!= nil{
		return err
	}
	limits[i], err = ConvertFromDefaultQuantity(unitSystem  , re.QUANTITIES[i], re.LIMITS[i]);
	if err!= nil{
		return err
	}
}
return RangeError{re.QUANTITIES, values, limits, unitSystem}
}	 


func (err RangeError)  Error() string{

	out := "";

   for i := 0; i < len(err.QUANTITIES); i++ {
			quantity :=err.QUANTITIES[i].String()
			   unit := units.GetUnit(err.QUANTITIES[i]);
			   message := "higher"
			   if err.VALUES[i] > err.LIMITS[i]{
				   message = "lower"
			   }
	   switch (i) {
		   case 0:
			   
			   quantity = strings.ToUpper(quantity[:1])+ quantity[:1];

			   out += fmt.Sprintf("%s value %g %s should be %s than %g %s",
			   quantity, err.VALUES[i], unit, message, err.LIMITS[i], unit);


		   default:
			   out += fmt.Sprintf(", when %s value %g %s is %s than %g %s",
					   quantity, err.VALUES[i], unit, message, err.LIMITS[i], unit);
	   }
   }
   out += ".";

   return out;
}		


// 	OutOfRangeException(IF97.Quantity[] quantities, double[] values, double[] limits, IF97.UnitSystem unitSystem) {

// 		if (quantities == null || values == null || limits == null) {
// 			throw new IllegalArgumentException("Arguments shouldn't be null.");

// 		} else if (quantities.length == 0 || values.length == 0 || limits.length == 0) {
// 			throw new IllegalArgumentException("Argument arrays shouldn't be empty.");

// 		} else if (quantities.length != values.length || values.length != limits.length) {
// 			throw new IllegalArgumentException("Argument arrays should have equal lengths.");
// 		}
// 		QUANTITIES = quantities.clone();
// 		VALUES = values.clone();
// 		LIMITS = limits.clone();
// 		UNIT_SYSTEM = unitSystem;
// 	}


