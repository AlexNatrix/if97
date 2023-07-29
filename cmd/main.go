package main

import (
	"fmt"

	"if97.com/cmd/lib/IF97"
)

func main() {
	var if97 = IF97.New(1)
	p,err:=if97.SaturationTemperatureP(20) //MPa cuz unit system is 1 -> IF97.New(1), T in Kelvins
	fmt.Print(p,err)
}
