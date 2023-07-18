package firstRegion_test

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"if97.com/cmd/lib/firstRegion"
)

func TestFirstRegion(t *testing.T) {
	cases := []struct {
		name   string
		f      func(pressure float64, temperature float64) float64
		values [][]float64
		tol    float64
	}{{
		name: "testHeatCapacityRatioPT",
		f:    firstRegion.HeatCapacityRatioPT,
		values: [][]float64{
			{0.417301218e1 / 0.412120160e1, 3, 300},
			{0.401008987e1 / 0.391736606e1, 80, 300},
			{0.465580682e1 / 0.322139223e1, 3, 500}},
		tol: 1e-8,
	}, {
		name: "testIsentropicExponentPT",
		f:    firstRegion.IsentropicExponentPT,
		values: [][]float64{
			{756.132, 3, 300},
			{34.212, 80, 298.15}, // 800 bar, 25 C, Table 3, p.282
			{34.394, 80, 300},
			{426.743, 3, 500},
		},
		tol: 1e-3,
	},

		{name: "testIsobaricCubicExpansionCoefficientPT",
			f: firstRegion.IsobaricCubicExpansionCoefficientPT,
			values: [][]float64{
				{0.277354533e-3, 3, 300},
				{0.344095843e-3, 80, 300},
				{0.164118128e-2, 3, 500}},
			tol: 1e-12,
		},
		{
			name: "testIsothermalCompressibilityPT",
			f:    firstRegion.IsothermalCompressibilityPT,
			values: [][]float64{
				{0.446382123e-3, 3, 300},
				{0.372039437e-3, 80, 300},
				{0.112892188e-2, 3, 500}},
			tol: 1e-11,
		},
		{
			name: "testSpecificEnthalpyPT",
			f:    firstRegion.SpecificEnthalpyPT,
			values: [][]float64{
				{0.115331273e3, 3, 300},
				{0.184142828e3, 80, 300},
				{0.975542239e3, 3, 500}},
			tol: 1e-6,
		},
		{
			name: "testSpecificEntropyPT",
			f:    firstRegion.SpecificEntropyPT,
			values: [][]float64{
				{0.392294792, 3, 300},
				{0.368563852, 80, 300},
				{0.258041912e1, 3, 500}},
			tol: 1e-9,
		},
		{
			name: "testSpecificInternalEnergyPT",
			f:    firstRegion.SpecificInternalEnergyPT,
			values: [][]float64{
				{0.112324818e3, 3, 300},
				{0.106448356e3, 80, 300},
				{0.971934985e3, 3, 500}},
			tol: 1e-6,
		},
		{
			name: "testSpecificIsobaricHeatCapacityPT",
			f:    firstRegion.SpecificIsobaricHeatCapacityPT,
			values: [][]float64{
				{0.417301218e1, 3, 300},
				{0.401008987e1, 80, 300},
				{0.465580682e1, 3, 500}},
			tol: 1e-8,
		},
		{
			name: "testSpecificIsochoricHeatCapacityPT",
			f:    firstRegion.SpecificIsochoricHeatCapacityPT,
			values: [][]float64{
				{0.412120160e1, 3, 300},
				{0.391736606e1, 80, 300},
				{0.322139223e1, 3, 500}},
			tol: 1e-8,
		},
		{
			name: "testSpecificVolumePT",
			f:    firstRegion.SpecificVolumePT,
			values: [][]float64{
				{0.100215168e-2, 3, 300},
				{0.971180894e-3, 80, 300},
				{0.120241800e-2, 3, 500}},
			tol: 1e-11,
		},

		{
			name: "testPressureHS",
			values: [][]float64{
				{9.800980612e-4, 0.001, 0}, // ok with tol = 1e-13
				{9.192954727e1, 90, 0},
				{5.868294423e1, 1500, 3.4}},
			f:   firstRegion.PressureHS,
			tol: 1e-8,
		}, {
			name: "testSpeedOfSoundPT",
			values: [][]float64{
				{0.150773921e4, 3, 300},
				{0.163469054e4, 80, 300},
				{0.124071337e4, 3, 500}},
			f:   firstRegion.SpeedOfSoundPT,
			tol: 1e-5,
		},
		{
			name: "testTemperaturePH",
			values: [][]float64{
				{0.391798509e3, 3, 500},
				{0.378108626e3, 80, 500},
				{0.611041229e3, 80, 1500}},
			f:   firstRegion.TemperaturePH,
			tol: 1e-6,
		},
		{
			name: "testTemperaturePS",
			values: [][]float64{
				{0.307842258e3, 3, 0.5},
				{0.309979785e3, 80, 0.5},
				{0.565899909e3, 80, 3}},
			f:   firstRegion.TemperaturePS,
			tol: 1e-6,
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
