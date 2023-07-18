package constants

const (
	R = 0.461526 //Specific gas constant of ordinary water [kJ/kg-K].

	Rm = 8.31451 //Molar gas constant of ordinary water [kJ/kmol-K].

	M = 18.015257 //Molar mass of ordinary water [kg/kmol].

	T0 = 273.15 //in most literature (it states 273.16 somewhere)

	Pc = 22.064 //Critical pressure [MPa].

	Tc = 647.096 //Critical temperature [K].

	Hc = 2087.546845 //Critical enthalpy [kJ/kg].

	Rhoc = 322 //Critical density [kg/m3].

	BTU = 1.055056 //British thermal unit acc. International standard ISO 31-4 on Quantities and units—Part 4: Heat, Appendix A [kJ]

	G = 9.80665 //Gravitational accelleration [m/s^2]

	Lb = 0.45359237 //International avoirdupois pound-mass [kg]

	Ft = 0.3048 // foot [m]

	FtSquare = Ft * Ft // square foot [m^2]

	FtCube = Ft * FtSquare // cubic foot [m^3]

	Hr = 3600 // hour [s]

	In = Ft / 12 // inch [m]

	InSquare = In * In // square inch [m^2]

	Lbf = Lb * G
	Psi = 1e-6 * Lbf / InSquare // pounds per square inch [MPa]

	Ra  = 5.0 / 9.0 // Rankine [K]
	T13 = T0 + 350
	//Critical entropy [kJ/kg-K].
)
