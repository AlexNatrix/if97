package thirdRegion

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestThirdRegion(t *testing.T) {
	cases := []struct {
		name   string
		f      func(pressure float64, temperature float64) float64
		values [][]float64
		tol    float64
	}{
		{name: "testHeatCapacityRatioRhoT",
			values: [][]float64{{0.138935717e2 / 0.319131787e1, 500, 650},
				{0.446579342e2 / 0.404118076e1, 200, 650},
				{0.634165359e1 / 0.271701677e1, 500, 750}},
			f:   REGION3.HeatCapacityRatioRhoT,
			tol: 1e-7,
		},
		{name: "testIsobaricCubicExpansionCoefficientRhoT",
			values: [][]float64{{0.168653107e-1, 500, 650},
				{0.685312229e-1, 200, 650},
				{0.441515098e-2, 500, 750}},
			f:   REGION3.IsobaricCubicExpansionCoefficientRhoT,
			tol: 1e-7,
		},
		{name: "testIsothermalCompressibilityRhoT",
			values: [][]float64{{0.345506956e-1, 500, 650},
				{0.375798565, 200, 650},
				{0.806710817e-2, 500, 750}},
			f:   REGION3.IsothermalCompressibilityRhoT,
			tol: 1e-9,
		},
		{name: "testIsothermalStressCoefficientRhoT",
			values: [][]float64{{0.565652647e3, 500, 650},
				{0.238728962e2, 200, 650},
				{0.791475213e3, 500, 750}},
			f:   REGION3.IsothermalStressCoefficientRhoT,
			tol: 1e-6,
		},
		{name: "testPressureHS",
			values: [][]float64{
				{2.555703246e1, 1700, 3.8},
				{4.540873468e1, 2000, 4.2},
				{6.078123340e1, 2100, 4.3},
				{6.363924887e1, 2400, 4.7},
				{3.434999263e1, 2600, 5.1},
				{8.839043281e1, 2700, 5.0}},
			f:   REGION3.PressureHS,
			tol: 1e-8,
		},
		{name: "testPressureRhoT",
			values: [][]float64{
				{0.255837018e2, 500, 650},
				{0.222930643e2, 200, 650},
				{0.783095639e2, 500, 750},
			},
			f:   REGION3.PressureRhoT,
			tol: 1e-7,
		},
		{name: "testSpecificEnthalpyRhoT",
			values: [][]float64{
				{0.186343019e4, 500, 650},
				{0.237512401e4, 200, 650},
				{0.225868845e4, 500, 750},
			},
			f:   REGION3.SpecificEnthalpyRhoT,
			tol: 1e-5,
		},
		{name: "testSpecificEntropyRhoT",
			values: [][]float64{
				{0.405427273e1, 500, 650},
				{0.485438792e1, 200, 650},
				{0.446971906e1, 500, 750},
			},
			f:   REGION3.SpecificEntropyRhoT,
			tol: 1e-8,
		},
		{name: "testSpecificInternalEnergyRhoT",
			values: [][]float64{
				{0.181226279e4, 500, 650},
				{0.226365868e4, 200, 650},
				{0.210206932e4, 500, 750},
			},
			f:   REGION3.SpecificInternalEnergyRhoT,
			tol: 1e-5,
		},
		{name: "testSpecificIsobaricHeatCapacityRhoT",
			values: [][]float64{
				{0.138935717e2, 500, 650},
				{0.446579342e2, 200, 650},
				{0.634165359e1, 500, 750},
			},
			f:   REGION3.SpecificIsobaricHeatCapacityRhoT,
			tol: 1e-7,
		},
		{name: "testSpecificIsochoricHeatCapacityRhoT",
			values: [][]float64{
				{0.319131787e1, 500, 650},
				{0.404118076e1, 200, 650},
				{0.271701677e1, 500, 750},
			},
			f:   REGION3.SpecificIsochoricHeatCapacityRhoT,
			tol: 1e-8,
		},
		{name: "testSpecificVolumePH",
			values: [][]float64{
				{1.749903962e-3, 20, 1700},
				{1.908139035e-3, 50, 2000},
				{1.676229776e-3, 100, 2100},
				{6.670547043e-3, 20, 2500},
				{2.801244590e-3, 50, 2400},
				{2.404234998e-3, 100, 2700},
			},
			f:   REGION3.SpecificVolumePH,
			tol: 1e-12,
		},
		{name: "testSpecificVolumePS",
			values: [][]float64{
				{1.733791463e-3, 20, 3.8},
				{1.469680170e-3, 50, 3.6},
				{1.555893131e-3, 100, 4.0},
				{6.262101987e-3, 20, 5.0},
				{2.332634294e-3, 50, 4.5},
				{2.449610757e-3, 100, 5.0},
			},
			f:   REGION3.SpecificVolumePS,
			tol: 1e-12,
		},
		{name: "testSpecificVolumePT",
			values: [][]float64{
				{1.470853100e-3, 50, 630}, // a
				{1.503831359e-3, 80, 670},
				{2.204728587e-3, 50, 710}, // b
				{1.973692940e-3, 80, 750},
				{1.761696406e-3, 20, 630}, // c
				{1.819560617e-3, 30, 650},
				{2.245587720e-3, 26, 656}, // d
				{2.506897702e-3, 30, 670},
				{2.970225962e-3, 26, 661}, // e
				{3.004627086e-3, 30, 675},
				{5.019029401e-3, 26, 671}, // f
				{4.656470142e-3, 30, 690},
				{2.163198378e-3, 23.6, 649}, // g
				{2.166044161e-3, 24, 650},
				{2.651081407e-3, 23.6, 652}, // h
				{2.967802335e-3, 24, 654},
				{3.273916816e-3, 23.6, 653}, // i
				{3.550329864e-3, 24, 655},
				{4.545001142e-3, 23.5, 655}, // j
				{5.100267704e-3, 24, 660},
				{6.109525997e-3, 23, 660}, // k
				{6.427325645e-3, 24, 670},
				{2.117860851e-3, 22.6, 646}, // l
				{2.062374674e-3, 23, 646},
				{2.533063780e-3, 22.6, 648.6}, // m
				{2.572971781e-3, 22.8, 649.3},
				{2.923432711e-3, 22.6, 649.0}, // n
				{2.913311494e-3, 22.8, 649.7},
				{3.131208996e-3, 22.6, 649.1}, // o
				{3.221160278e-3, 22.8, 649.9},
				{3.715596186e-3, 22.6, 649.4}, // p
				{3.664754790e-3, 22.8, 650.2},
				{1.970999272e-3, 21.1, 640}, // q
				{2.043919161e-3, 21.8, 643},
				{5.251009921e-3, 21.1, 644}, // r
				{5.256844741e-3, 21.8, 648},
				{1.932829079e-3, 19.1, 635}, // s
				{1.985387227e-3, 20, 638},
				{8.483262001e-3, 17, 626}, // t
				{6.227528101e-3, 20, 640},
				{2.268366647e-3, 21.5, 644.6}, // u
				{2.296350553e-3, 22, 646.1},
				{2.832373260e-3, 22.5, 648.6}, // v
				{2.811424405e-3, 22.3, 647.9},
				{3.694032281e-3, 22.15, 647.5}, // w
				{3.622226305e-3, 22.3, 648.1},
				{4.528072649e-3, 22.11, 648}, // x
				{4.556905799e-3, 22.3, 649},
				{2.698354719e-3, 22, 646.84}, // y
				{2.717655648e-3, 22.064, 647.05},
				{3.798732962e-3, 22, 646.89}, // z
				{3.701940010e-3, 22.064, 647.15},
			},
			f:   REGION3.SpecificVolumePT,
			tol: 1e-12,
		},
		{name: "testSpeedOfSoundRhoT",
			values: [][]float64{
				{0.502005554e3, 500, 650},
				{0.383444594e3, 200, 650},
				{0.760696041e3, 500, 750},
			},
			f:   REGION3.SpeedOfSoundRhoT,
			tol: 1e-6,
		},
		{name: "testTemperaturePH",
			values: [][]float64{
				{6.293083892e2, 20, 1700},
				{6.905718338e2, 50, 2000},
				{7.336163014e2, 100, 2100},
				{6.418418053e2, 20, 2500},
				{7.351848618e2, 50, 2400},
				{8.420460876e2, 100, 2700},
			},
			f:   REGION3.TemperaturePH,
			tol: 1e-7,
		},
		{name: "testTemperaturePS",
			values: [][]float64{
				{6.282959869e2, 20, 3.8},
				{6.297158726e2, 50, 3.6},
				{7.056880237e2, 100, 4.0},
				{6.401176443e2, 20, 5.0},
				{7.163687517e2, 50, 4.5},
				{8.474332825e2, 100, 5.0},
			},
			f:   REGION3.TemperaturePS,
			tol: 1e-7,
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
	t.Run("testIsentropicExponentRhoT", func(t *testing.T) {
		for _, x := range [][]float64{{6.5889, 50, 673.15},
			{8.0282, 80, 673.15},
			{3.0622, 80, 773.15}} {
			rho := 1 / REGION3.SpecificVolumePT(x[1], x[2])
			assert.InDelta(t, x[0], REGION3.IsentropicExponentRhoT(rho, x[2]), 1e-4)
		}
	})
	t.Run("testEnthalpy2bc", func(t *testing.T) {
		for _, x := range [][]float64{{2.095936454e3, 25}} {
			assert.InDelta(t, x[0], enthalpy3ab(x[1]), 1e-6)
		}
	})
}
