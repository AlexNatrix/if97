package IF97

type UNITSYSTEM int

const (
	_ UNITSYSTEM = iota
	DEFAULT
)

type WSPcalculator interface {
}

type IF97 struct {
	UnitSystem int
}

func convertToDefault(quantity []float64, value float64) float64 {
	return value*quantity[0] + quantity[1]
}


func PrandtlHS(enthalpy float64, entropy float64) (float64 , error) {

	h := convertToDefault(UNIT_SYSTEM.SPECIFIC_ENTHALPY, enthalpy)
	s := convertToDefault(UNIT_SYSTEM.SPECIFIC_ENTROPY, entropy)

	
	region,err = Region.getRegionHS(h, s);
	if err != nil{
		return -1,err
	}
	p := region.pressureHS(h, s)
	T := region.temperatureHS(h, s)

	return Calculate.PrandtlPT(p, T),nil;
}


func PrandtlPH(double pressure, double enthalpy) (float64, error) {

	 p := convertToDefault(UNIT_SYSTEM.PRESSURE, pressure)
	 h := convertToDefault(UNIT_SYSTEM.SPECIFIC_ENTHALPY, enthalpy)

	try {
		return calculator.PrandtlPH(p, h);

	} catch (OutOfRangeException e) {
		throw e.convertFromDefault(UNIT_SYSTEM);
	}
}
