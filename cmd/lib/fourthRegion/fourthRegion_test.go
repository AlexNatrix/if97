package fourthRegion_test

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"if97.com/cmd/lib/fourthRegion"
	"if97.com/cmd/lib/utils/constants"
)

func TestFourthRegion(t *testing.T) {

	t.Run("testSpecificVolumeSaturatedLiquidT", func(t *testing.T) {
		t.Parallel()
		for _, x := range [][]float64{
			{0.00229020, 644.15},
			{0.00238170, 645.15},
			{0.00252643, 646.15}} {
			ps := fourthRegion.SaturationPressureT(x[1])
			assert.InDelta(t, x[0], fourthRegion.SpecificVolumeSaturatedLiquidP(ps), 1e-8)
		}
	})

	t.Run("testSpecificEnthalpyTX", func(t *testing.T) {
		t.Parallel()
		for _, x := range [][]float64{
			{1686.2747, 625, 0},
			{1798.6443, 625, 0.13},
			{2118.4653, 625, 0.5},
			{2542.0121, 625, 0.99},
			{2550.6559, 625, 1}} {
			ps := fourthRegion.SaturationPressureT(x[1])
			assert.InDelta(t, x[0], fourthRegion.SpecificEnthalpyPX(ps, x[2]), 1e-2)
		}
	})

	casesArityOne := []struct {
		name   string
		f      func(pressure float64) float64
		values [][]float64
		tol    float64
	}{
		{
			name: "testSaturationPressureB34H",
			f:    fourthRegion.SaturationPressureB34H,
			values: [][]float64{
				{1.724175718e1, 1700},
				{2.193442957e1, 2000},
				{2.018090839e1, 2400},
			},
			tol: 1e-8,
		},
		{
			name: "testSaturationPressureB34S",
			f:    fourthRegion.SaturationPressureB34S,
			values: [][]float64{
				{1.687755057e1, 3.8},
				{2.164451789e1, 4.2},
				{1.668968482e1, 5.2},
			},
			tol: 1e-8,
		},
		{
			name: "testSaturationPressureT",
			f:    fourthRegion.SaturationPressureT,
			values: [][]float64{
				{0.353658941e-2, 300},
				{0.263889776e1, 500},
				{0.123443146e2, 600},
			},
			tol: 1e-7,
		},
		{
			name: "testSaturationTemperatureP",
			f:    fourthRegion.SaturationTemperatureP,
			values: [][]float64{
				{0.372755919e3, 0.1},
				{0.453035632e3, 1},
				{0.584149488e3, 10},
			},
			tol: 1e-6,
		},
		{
			name: "testSpecificEnthalpySaturatedLiquidP",
			f:    fourthRegion.SpecificEnthalpySaturatedLiquidP,
			values: [][]float64{
				{1407.87, 10},
				{1610.15, 15},
				{1827.10, 20},
				{1855.90, 20.5},
				{1889.40, 21},
				{1932.81, 21.5},
				{2021.92, 22},
				//FAIL {1933.0015, 21.5}, // Zittau
				//FAIL {2013.3573, 22}, // Zittau
				{2087.55, constants.Pc},
			},
			tol: 1e-2,
		},
		{
			name: "testSpecificEnthalpySaturatedVapourP",
			f:    fourthRegion.SpecificEnthalpySaturatedVapourP,
			values: [][]float64{
				{2725.47, 10},
				{2610.86, 15},
				{2411.39, 20},
				{2378.16, 20.5},
				{2337.54, 21},
				{2282.18, 21.5},
				//FAIL {2164.18, 22},
				//FAIL {2281.8484, 21.5}, // Zittau
				//FAIL {2163.2117, 22}, // Zittau
				{2087.55, constants.Pc},
			},
			tol: 1e-2,
		},
		{
			name: "testSpecificEntropySaturatedLiquidP",
			f:    fourthRegion.SpecificEntropySaturatedLiquidP,
			values: [][]float64{
				{3.36029, 10},
				{3.68445, 15},
				{4.01538, 20},
				{4.0588, 20.5},
				{4.10926, 21},
				{4.1749, 21.5}, // Wagner Table 2
				{4.3109, 22},   // Wagner Table 2
				//FAIL {4.17520, 21.5}, // Zittau
				//FAIL {4.29763, 22}, // Zittau
				{4.41202, constants.Pc},
			},
			tol: 1e-4,
		},
		{
			name: "testSpecificEntropySaturatedVapourP",
			f:    fourthRegion.SpecificEntropySaturatedVapourP,
			values: [][]float64{
				{5.61589, 10},
				{5.31080, 15},
				{4.92991, 20},
				{4.87357, 20.5},
				{4.80624, 21},
				{4.7166, 21.5}, // Wagner Table 2
				{4.5308, 22},   // Wagner Table 2
				//FAIL {4.71608, 21.5}, // Zittau
				//FAIL {4.5293, 22}, // Zittau
				{4.41202, constants.Pc},
			},
			tol: 1e-4,
		},
		{
			name: "testSpecificVolumeSaturatedLiquidP",
			f:    fourthRegion.SpecificVolumeSaturatedLiquidP,
			values: [][]float64{
				{0.00145262, 10},
				{0.00165696, 15},
				{0.00168249, 15.5},
				{0.00170954, 16},
				{0.00173833, 16.5},
				{0.00176934, 17},
				{0.00203865, 20},
				{0.00211358, 20.5},
				{0.00221186, 21},
				{0.00236016, 21.5},
				{0.00275039, 22},
				//FAIL {0.00270572, 22}, // Zittau
				{0.00310559, constants.Pc},
			},
			tol: 1e-8,
		},
		{
			name: "testSpecificVolumeSaturatedVapourP",
			f:    fourthRegion.SpecificVolumeSaturatedVapourP,
			values: [][]float64{
				{0.0103401, 15},
				{0.00981114, 15.5},
				{0.00930813, 16},
				{0.00882826, 16.5},
				{0.00836934, 17},
				{0.00585828, 20},
				{0.00543778, 20.5},
				{0.00498768, 21},
				{0.00446300, 21.5},
				{0.00357662, 22},
				{0.00310559, constants.Pc},
			},
			tol: 1e-8,
		},
	}
	for _, tc := range casesArityOne {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			for _, x := range tc.values {
				assert.InDelta(t, x[0], tc.f(x[1]), tc.tol)
			}
		})
	}

	// //    @Test
	//     public void testSpecificVolumeSaturatedVapour_fail16_53() {

	//         double[][] X = {
	//             {0.00882826, 16.5},
	//             {0.00882826, 16.51},
	//             {0.00882826, 16.52},
	//             {0.00882826, 16.529},
	//             {0.00882826, 16.53}, //XXX Problematic value
	//             {0.00882826, 16.531},
	//             {0.00882826, 16.54},
	//             {0.00882826, 16.55}
	//         };

	//         for (double[] x : X) {
	//             double Tsat = region4.saturationTemperatureP(x[1]),
	//                     nu = region4.specificVolumeSaturatedVapourP(x[1]),
	//                     nu2 = Region.REGION2.specificVolumePT(x[1], Tsat),
	//                     nu3 = Region.REGION3.specificVolumePT(x[1], Tsat);
	//             System.out.format("p: %.3f, Ts: %.3f", x[1], Tsat);
	//             System.out.format(", specific volume: %8.6f", nu);
	//             System.out.format(", region 2: %8.6f", nu2);
	//             System.out.format(", region 3: %8.6f%n", nu3);
	//         }
	//     }

	cases := []struct {
		name   string
		f      func(pressure float64, temperature float64) float64
		values [][]float64
		tol    float64
	}{{name: " testSpecificEnthalpyPX",
		values: [][]float64{
			{1827.1005, 20, 0},
			{1903.0579, 20, 0.13},
			{2119.2443, 20, 0.5},
			{2405.5451, 20, 0.99},
			{2411.3880, 20, 1}},
		f:   fourthRegion.SpecificEnthalpyPX,
		tol: 1e-3,
	},
		{name: "testSpecificVolumePX",
			values: [][]float64{
				{0.00145262, 10, 0},
				{0.00203865, 20, 0},
				{0.0025352, 20, 0.13},
				{0.00394847, 20, 0.5},
				{0.00582009, 20, 0.99},
				{0.00585828, 20, 1},
				{0.0180336, 10, 1}},
			f:   fourthRegion.SpecificVolumePX,
			tol: 1e-7,
		},
		{name: "testVapourFractionPH",
			values: [][]float64{
				{0.00496239, 20, 1830},
				{0.295915, 20, 2000},
				{0.966917, 16, 2550},
				{0.980510, 20, 2400},
				{0.999641, 10, 2725}},
			f:   fourthRegion.VapourFractionPH,
			tol: 1e-6,
		},
		{name: "testVapourFractionPS",
			values: [][]float64{
				{0.681062, 5, 5},
				{0.726951, 10, 5},
				{0.420567, 20, 4.4},
				{0, 20, 4.01538},
				{0, 21, 4.10926},
				{0, 22, 4.31085},
				{1, 20, 4.92991},
				{1, 21, 4.80624},
				{1, 22, 4.5308}},
			f:   fourthRegion.VapourFractionPS,
			tol: 1e-4,
		},
		{name: "testTemperatureHS",
			values: [][]float64{
				//{Double.NaN, 1800, 5.2}, // should throw an out-of-range exception
				{3.468475498e2, 1800, 5.3},
				{4.251373305e2, 2400, 6},
				{5.225579013e2, 2500, 5.5}},
			f:   fourthRegion.TemperatureHS,
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

}
