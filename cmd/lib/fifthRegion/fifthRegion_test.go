package fifthRegion_test

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"if97.com/cmd/lib/fifthRegion"
)

func TestFifthRegion(t *testing.T) {
	t.Run("testIsentropicExponentPT", func(t *testing.T) {
		t.Parallel()
		for _, x := range [][]float64{
			{0.917068690e3, 0.5, 1500},
			{0.928548002e3, 30, 1500},
			{0.106736948e4, 30, 2000}} {
			w2 := x[0] * x[0] // speed of sound [m/s] squared [m2/s2]
			p := x[1] * 1e6   // [Pa]
			rho := 1 / fifthRegion.SpecificVolumePT(x[1], x[2])
			kappa := w2 * rho / p // [-]
			assert.InDelta(t, kappa, fifthRegion.IsentropicExponentPT(x[1], x[2]), 1e-8)
		}
	})

	cases := []struct {
		name   string
		f      func(pressure float64, temperature float64) float64
		values [][]float64
		tol    float64
	}{
		{
			name: "testHeatCapacityRatioPT",
			f:    fifthRegion.HeatCapacityRatioPT,
			values: [][]float64{
				{0.261609445e1 / 0.215337784e1, 0.5, 1500},
				{0.272724317e1 / 0.219274829e1, 30, 1500},
				{0.288569882e1 / 0.239589436e1, 30, 2000}},
			tol: 1e-8,
		},

		{
			name: "testIsobaricCubicExpansionCoefficientPT",
			f:    fifthRegion.IsobaricCubicExpansionCoefficientPT,
			values: [][]float64{
				{0.667539000e-3, 0.5, 1500},
				{0.716950754e-3, 30, 1500},
				{0.508830641e-3, 30, 2000}},
			tol: 1e-12,
		},
		{
			name: "testIsothermalCompressibilityPT",
			f:    fifthRegion.IsothermalCompressibilityPT,
			values: [][]float64{
				{0.200003859e1, 0.5, 1500},
				{0.332881253e-1, 30, 1500},
				{0.329193892e-1, 30, 2000}},
			tol: 1e-8,
		},
		{
			name: "testSpecificEnthalpyPT",
			f:    fifthRegion.SpecificEnthalpyPT,
			values: [][]float64{
				{0.521976855e4, 0.5, 1500},
				{0.516723514e4, 30, 1500},
				{0.657122604e4, 30, 2000}},
			tol: 1e-5,
		},
		{
			name: "testSpecificEntropyPT",
			f:    fifthRegion.SpecificEntropyPT,
			values: [][]float64{
				{0.965408875e1, 0.5, 1500},
				{0.772970133e1, 30, 1500},
				{0.853640523e1, 30, 2000}},
			tol: 1e-8,
		},
		{
			name: "testSpecificInternalEnergyPT",
			f:    fifthRegion.SpecificInternalEnergyPT,
			values: [][]float64{
				{0.452749310e4, 0.5, 1500},
				{0.447495124e4, 30, 1500},
				{0.563707038e4, 30, 2000}},
			tol: 1e-5,
		},
		{
			name: "testSpecificIsobaricHeatCapacityPT",
			f:    fifthRegion.SpecificIsobaricHeatCapacityPT,
			values: [][]float64{
				{0.261609445e1, 0.5, 1500},
				{0.272724317e1, 30, 1500},
				{0.288569882e1, 30, 2000}},
			tol: 1e-8,
		},
		{
			name: "testSpecificIsochoricHeatCapacityPT",
			f:    fifthRegion.SpecificIsochoricHeatCapacityPT,
			values: [][]float64{
				{0.215337784e1, 0.5, 1500},
				{0.219274829e1, 30, 1500},
				{0.239589436e1, 30, 2000}},
			tol: 1e-8,
		},
		{
			name: "testSpecificVolumePT",
			f:    fifthRegion.SpecificVolumePT,
			values: [][]float64{
				{0.138455090e1, 0.5, 1500},
				{0.230761299e-1, 30, 1500},
				{0.311385219e-1, 30, 2000}},
			tol: 1e-8,
		},
		{
			name: "testSpeedOfSoundPT",
			f:    fifthRegion.SpeedOfSoundPT,
			values: [][]float64{
				{0.917068690e3, 0.5, 1500},
				{0.928548002e3, 30, 1500},
				{0.106736948e4, 30, 2000}},
			tol: 1e-5,
		},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			for _, x := range tc.values {
				assert.InDelta(t, x[0], tc.f(x[1], x[2]), tc.tol)
			}
		})
	}
}
