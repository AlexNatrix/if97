package secondRegionMeta_test

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"if97.com/cmd/lib/secondRegionMeta"
)

func TestSecondRegion(t *testing.T) {
	cases := []struct {
		name   string
		f      func(pressure float64, temperature float64) float64
		values [][]float64
		tol    float64
	}{
		{
			name: "testIsobaricCubicExpansionCoefficientPT",
			values: [][]float64{
				{0.318819824e-2, 1, 450},
				{0.348506136e-2, 1, 440},
				{0.418276571e-2, 1.5, 450}},
			f:   secondRegionMeta.IsobaricCubicExpansionCoefficientPT,
			tol: 1e-11,
		},
		{
			name: "testIsothermalCompressibilityPT",
			values: [][]float64{
				{0.109364239e1, 1, 450},
				{0.111133230e1, 1, 440},
				{0.787967952, 1.5, 450}},
			f:   secondRegionMeta.IsothermalCompressibilityPT,
			tol: 1e-8,
		},
		{
			name: "testSpecificEnthalpyPT",
			values: [][]float64{
				{0.276881115e4, 1, 450},
				{0.274015123e4, 1, 440},
				{0.272134539e4, 1.5, 450}},
			f:   secondRegionMeta.SpecificEnthalpyPT,
			tol: 1e-5,
		},
		{
			name: "testSpecificEntropyPT",
			values: [][]float64{
				{0.656660377e1, 1, 450},
				{0.650218759e1, 1, 440},
				{0.629170440e1, 1.5, 450}},
			f:   secondRegionMeta.SpecificEntropyPT,
			tol: 1e-8,
		},
		{
			name: "testSpecificInternalEnergyPT",
			values: [][]float64{
				{0.257629461e4, 1, 450},
				{0.255393894e4, 1, 440},
				{0.253881758e4, 1.5, 450}},
			f:   secondRegionMeta.SpecificInternalEnergyPT,
			tol: 1e-5,
		},
		{
			name: "testSpecificIsobaricHeatCapacityPT",
			values: [][]float64{
				{0.276349265e1, 1, 450},
				{0.298166443e1, 1, 440},
				{0.362795578e1, 1.5, 450}},
			f:   secondRegionMeta.SpecificIsobaricHeatCapacityPT,
			tol: 1e-8,
		},
		{
			name: "testSpecificIsochoricHeatCapacityPT",
			values: [][]float64{
				{0.195830730e1, 1, 450},
				{0.208622142e1, 1, 440},
				{0.241213708e1, 1.5, 450}},
			f:   secondRegionMeta.SpecificIsochoricHeatCapacityPT,
			tol: 1e-8,
		},
		{
			name: "testSpecificVolumePT",
			values: [][]float64{
				{0.192516540, 1, 450},
				{0.186212297, 1, 440},
				{0.121685206, 1.5, 450}},
			f:   secondRegionMeta.SpecificVolumePT,
			tol: 1e-9,
		},
		{
			name: "testSpeedOfSoundPT",
			values: [][]float64{{0.498408101e3, 1, 450},
				{0.489363295e3, 1, 440},
				{0.481941819e3, 1.5, 450},
			},
			f:   secondRegionMeta.SpeedOfSoundPT,
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
