**IF97 - water and steam properties**

_Steam tables for industrial use according to the international standard for the properties of water and steam,
the IAPWS-IF97 formulation and the international standards for transport and other properties._

__Use:__
```
package main

import (
	"fmt"

	"if97.com/cmd/lib/IF97"
)

func main() {
	var if97 = IF97.New(1)
	p,err:=if97.SaturationTemperatureP(20) //MPa cuz unit system is 1 -> IF97.New(1), T in Kelvins
	fmt.Print(p,err) //638.8959115457051 <nil>
}
```

_Default unit system is_
```
{
    "MPa",
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
```

__Read docs for more explanation.__
