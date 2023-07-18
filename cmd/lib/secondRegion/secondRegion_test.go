package secondRegion_test

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"if97.com/cmd/lib/secondRegion"
)

func TestSecondRegion(t *testing.T) {
	cases := []struct {
		name   string
		f      func(pressure float64, temperature float64) float64
		values [][]float64
		tol    float64
	}{{name: "testHeatCapacityRatioPT",
		values: [][]float64{
			{0.191300162e1 / 0.144132662e1, 0.0035, 300},
			{0.208141274e1 / 0.161978333e1, 0.0035, 700},
			{0.103505092e2 / 0.297553837e1, 30, 700}},
		f:   secondRegion.HeatCapacityRatioPT,
		tol: 1e-8,
	},
		{
			name: "testIsentropicExponentPT",
			values: [][]float64{
				{1.2881, 0.1, 673.15},
				{1.2935, 20, 673.15},
				{1.4227, 50, 873.15}},
			f:   secondRegion.IsentropicExponentPT,
			tol: 1e-4,
		},
		{
			name: "testIsobaricCubicExpansionCoefficientPT",
			values: [][]float64{
				{0.337578289e-2, 0.0035, 300},
				{0.142878736e-2, 0.0035, 700},
				{0.126019688e-1, 30, 700}},
			f:   secondRegion.IsobaricCubicExpansionCoefficientPT,
			tol: 1e-10,
		},
		{
			name: "testIsothermalCompressibilityPT",
			values: [][]float64{
				{0.286239651e3, 0.0035, 300},
				{0.285725461e3, 0.0035, 700},
				{0.818411389e-1, 30, 700}},
			f:   secondRegion.IsothermalCompressibilityPT,
			tol: 1e-6,
		},
		{
			name: "testPressureHS",
			values: [][]float64{
				{1.371012767, 2800, 6.5},
				{1.879743844e-3, 2800, 9.5},
				{1.024788997e-1, 4100, 9.5},
				{4.793911442, 2800, 6},
				{8.395519209e1, 3600, 6},
				{7.527161441, 3600, 7},
				{9.439202060e1, 2800, 5.1},
				{8.414574124, 2800, 5.8},
				{8.376903879e1, 3400, 5.8}},
			f:   secondRegion.PressureHS,
			tol: 1e-8,
		},
		{
			name: "testSpecificEnthalpyPT",
			values: [][]float64{
				{0.254991145e4, 0.0035, 300},
				{0.333568375e4, 0.0035, 700},
				{0.263149474e4, 30, 700}},
			f:   secondRegion.SpecificEnthalpyPT,
			tol: 1e-5,
		},
		{
			name: "testSpecificEntropyPT",
			values: [][]float64{
				{0.852238967e1, 0.0035, 300},
				{0.101749996e2, 0.0035, 700},
				{0.517540298e1, 30, 700}},
			f:   secondRegion.SpecificEntropyPT,
			tol: 1e-7,
		},
		{
			name: "testSpecificInternalEnergyPT",
			values: [][]float64{
				{0.241169160e4, 0.0035, 300},
				{0.301262819e4, 0.0035, 700},
				{0.246861076e4, 30, 700}},
			f:   secondRegion.SpecificInternalEnergyPT,
			tol: 1e-5,
		},
		{
			name: "testSpecificIsobaricHeatCapacityPT",
			values: [][]float64{
				{0.191300162e1, 0.0035, 300},
				{0.208141274e1, 0.0035, 700},
				{0.103505092e2, 30, 700}},
			f:   secondRegion.SpecificIsobaricHeatCapacityPT,
			tol: 1e-8,
		},
		{
			name: "testSpecificIsochoricHeatCapacityPT",
			values: [][]float64{
				{0.144132662e1, 0.0035, 300},
				{0.161978333e1, 0.0035, 700},
				{0.297553837e1, 30, 700}},
			f:   secondRegion.SpecificIsochoricHeatCapacityPT,
			tol: 1e-8,
		},
		{
			name: "testSpecificVolumePT",
			values: [][]float64{
				{0.394913866e2, 0.0035, 300},
				{0.923015898e2, 0.0035, 700},
				{0.542946619e-2, 30, 700}},
			f:   secondRegion.SpecificVolumePT,
			tol: 1e-4,
		},
		{
			name: "testSpeedOfSoundPT",
			values: [][]float64{
				{0.427920172e3, 0.0035, 300},
				{0.644289068e3, 0.0035, 700},
				{0.480386523e3, 30, 700}},
			f:   secondRegion.SpeedOfSoundPT,
			tol: 1e-6,
		},
		{
			name: "testTemperaturePH",
			values: [][]float64{
				{0.534433241e3, 0.001, 3000},
				{0.575373370e3, 3, 3000},
				{0.101077577e4, 3, 4000},
				{0.801299102e3, 5, 3500},
				{0.101531583e4, 5, 4000},
				{0.875279054e3, 25, 3500},
				{0.743056411e3, 40, 2700},
				{0.791137067e3, 60, 2700},
				{0.882756860e3, 60, 3200}},
			f:   secondRegion.TemperaturePH,
			tol: 1e-5,
		},
		{
			name: "testTemperaturePS",
			values: [][]float64{
				{0.399517097e3, 0.1, 7.5},
				{0.514127081e3, 0.1, 8},
				{0.103984917e4, 2.5, 8},
				{0.600484040e3, 8, 6},
				{0.106495556e4, 8, 7.5},
				{0.103801126e4, 90, 6},
				{0.697992849e3, 20, 5.75},
				{0.854011484e3, 80, 5.25},
				{0.949017998e3, 80, 5.75}},
			f:   secondRegion.TemperaturePS,
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
