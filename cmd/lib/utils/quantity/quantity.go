package quantity

const (

// = Quantity.getPartialDerivatives();
)

// Quantities as defined by reference given above.
type Quantity struct {
	NAME   string
	SYMBOL string
}

var (
	P = Quantity{"absolute pressure", "p"} //Absolute pressure.

	T = Quantity{"temperature", "T"} // Temperature

	V = Quantity{"specific volume", "v"} //Specific volume.

	U = Quantity{"specific internal energy", "u"} //Specific internal energy.

	H = Quantity{"specific enthalpy", "h"} //Specific enthalpy.

	S = Quantity{"specific entropy", "s"} //Specific entropy.

	G = Quantity{"specific Gibbs free energy", "g"} //Specific Gibbs free energy.

	F = Quantity{"specific Helmholtz free energy", "f"} //Specific Helmholtz free energy.

	Rho = Quantity{"density", "\u03c1"} //Density.

	A = Quantity{"thermal diffusivity", "a"} //Thermal diffusivity.

	Cp = Quantity{"specific isobaric heat capacity", "<html>c<sub>p</sub></html>"} //Specific isobaric heat capacity.

	Cv = Quantity{"specific isochoric heat capacity", "<html>c<sub>v</sub></html>"} //Specific isochoric heat capacity.

	N = Quantity{"refractive index", "n"} //Refractive index.

	W = Quantity{"speed of sound", "w"} //Speed of sound.

	X = Quantity{"vapour fraction", "x"} //Vapour fraction.

	AlphaV = Quantity{"isobaric cubic expansion coefficient", "<html>&alpha;<sub>v</sub></html>"} //Isobaric cubic expansion coefficient.

	Epsilon = Quantity{"dielectric constant", "\u03b5"} //Dielectric constant.

	Eta = Quantity{"dynamic viscosity", "\u03bc"} //Dynamic viscosity.

	Kappa = Quantity{"thermal diffusivity", "\u03ba"} //Thermal diffusivity.

	KappaT = Quantity{"isothermal compressibility", "<html>&kappa;<sub>T</sub></html>"} //Isothermal compressibility

	Lambda = Quantity{"thermal conductivity", "\u03bb"} //Thermal conductivity.

	LambdaL = Quantity{"wavelength", "<html>&lambda;<sub>L</sub></html>"} //Wavelength of light.

	Nu = Quantity{"kinematic viscosity", "\u03bd"} //Kinematic viscosity.

	Sigma = Quantity{"surface tension", "\u03c3"} //Surface tension.

	Pr = Quantity{"Prandtl number", "Pr"} //Prandtl number.

	Z = Quantity{"compression factor", "Z"} //Compression factor or real-gas factor.

	Gamma = Quantity{"isentropic exponent", "\u03b3"} //Isentropic exponent, or heat capacity ratio.

)

func getPartialDerivatives() []Quantity {
	return []Quantity{
		P, T, V, U, H, S, G, F, Rho,
	}
}

func quantitySet() []Quantity {
	return []Quantity{
		P, T, V, U, H, S, G, F, Rho, A, Cp, Cv, N,
		W,
		X,
		AlphaV,
		Epsilon,
		Eta,
		Kappa,
		KappaT,
		Lambda,
		LambdaL,
		Nu,
		Sigma,
		Pr,
		Z,
		Gamma,
	}
}

func (q *Quantity) String() string {
	return q.NAME
}
