package IF97

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"if97.com/cmd/lib/thirdRegion"
)

func TestIF97(t *testing.T) {
	var if97 = New(1)

	cases := []struct {
		f      func(p float64, q float64) (float64, error)
		name   string
		values [][]float64
		tol    float64
	}{
		{
			name: "HeatCapacityRatioPT",
			f:    if97.HeatCapacityRatioPT,
			values: [][]float64{
				{0.417301218e1 / 0.412120160e1, 3, 300}, // region 1
				{0.401008987e1 / 0.391736606e1, 80, 300},
				{0.465580682e1 / 0.322139223e1, 3, 500},
				{0.191300162e1 / 0.144132662e1, 0.0035, 300}, // region 2
				{0.208141274e1 / 0.161978333e1, 0.0035, 700},
				{0.103505092e2 / 0.297553837e1, 30, 700},
				{0.261609445e1 / 0.215337784e1, 0.5, 1500}, // region 5
				{0.272724317e1 / 0.219274829e1, 30, 1500},
				{0.288569882e1 / 0.239589436e1, 30, 2000}},
			tol: 1e-8,
		},
		{
			name: "testIsentropicExponentPT",
			f:    if97.IsentropicExponentPT,
			values: [][]float64{
				{756.132, 3, 300},    // region 1
				{34.212, 80, 298.15}, // 800 bar, 25 C, Table 3, p.282
				{34.394, 80, 300},
				{426.743, 3, 500},
				{1.2881, 0.1, 673.15}, // region 2
				{1.2935, 20, 673.15},
				{1.4227, 50, 873.15}},
			tol: 1e-3,
		},
		{
			name: "testPressureHS",
			f:    if97.PressureHS,
			values: [][]float64{
				{9.800980612e-4, 0.001, 0}, // region 1: ok with tol = 1e-13
				{9.192954727e1, 90, 0},
				{5.868294423e1, 1500, 3.4},
				{1.371012767, 2800, 6.5}, // region 2
				{1.879743844e-3, 2800, 9.5},
				{1.024788997e-1, 4100, 9.5},
				{4.793911442, 2800, 6},
				{8.395519209e1, 3600, 6},
				{7.527161441, 3600, 7},
				{9.439202060e1, 2800, 5.1},
				{8.414574124, 2800, 5.8},
				{8.376903879e1, 3400, 5.8},
				{2.555703246e1, 1700, 3.8}, // region 3
				{4.540873468e1, 2000, 4.2},
				{6.078123340e1, 2100, 4.3},
				{6.363924887e1, 2400, 4.7},
				{3.434999263e1, 2600, 5.1},
				{8.839043281e1, 2700, 5.0}},
			tol: 1e-8,
		},
		{
			name: "testSpecificEnthalpyPT",
			f:    if97.SpecificEnthalpyPT,
			values: [][]float64{
				{0.115331273e3, 3, 300}, // region 1
				{0.184142828e3, 80, 300},
				{0.975542239e3, 3, 500},
				{0.254991145e4, 0.0035, 300}, // region 2 ==> actually metastable region
				{0.333568375e4, 0.0035, 700}, // region 2 ==> actually metastable region
				{0.263149474e4, 30, 700},     // region 2
				//{0.276881115e4, 1, 450}, // region 2 metastable-vapour ==> actually region 1
				//{0.274015123e4, 1, 440},
				//{0.272134539e4, 1.5, 450},
				{0.521976855e4, 0.5, 1500}, // region 5
				{0.516723514e4, 30, 1500},
				{0.657122604e4, 30, 2000}},
			tol: 1e-5,
		},
		/**
		 * Tests specific entropy as a function of pressure and temperature.
		 *
		 * Disabled region 2 and region 2 meta tests aren't actually in these
		 * regions.
		 */
		{
			name: "testSpecificEntropyPT",
			f:    if97.SpecificEntropyPT,
			values: [][]float64{
				{0.392294792, 3, 300}, // region 1 tests
				{0.368563852, 80, 300},
				{0.258041912e1, 3, 500},
				//{0.852238967e1, 0.0035, 300}, // region 2 tests
				//{0.101749996e2, 0.0035, 700},
				{0.517540298e1, 30, 700},
				//{0.656660377e1, 1, 450}, // region 2 meta tests
				//{0.650218759e1, 1, 440},
				//{0.629170440e1, 1.5, 450},
				{0.965408875e1, 0.5, 1500}, // region 5 tests
				{0.772970133e1, 30, 1500},
				{0.853640523e1, 30, 2000}},
			tol: 1e-5,
		},

		{
			name: "testSpecificVolumePH",
			f:    if97.SpecificVolumePH,
			values: [][]float64{
				{1.749903962e-3, 20, 1700}, // region 3
				{1.908139035e-3, 50, 2000},
				{1.676229776e-3, 100, 2100},
				{6.670547043e-3, 20, 2500},
				{2.801244590e-3, 50, 2400},
				{2.404234998e-3, 100, 2700}},
			tol: 1e-12,
		},

		{
			name: "testSpecificVolumePS",
			f:    if97.SpecificVolumePS,
			values: [][]float64{
				{1.733791463e-3, 20, 3.8}, // region 3
				{1.469680170e-3, 50, 3.6},
				{1.555893131e-3, 100, 4.0},
				{6.262101987e-3, 20, 5.0},
				{2.332634294e-3, 50, 4.5},
				{2.449610757e-3, 100, 5.0}},
			tol: 1e-12,
		},

		{
			name: "testSpecificVolumePT",
			f:    if97.SpecificVolumePT,
			values: [][]float64{
				{0.100215168e-2, 3, 300}, // region 1
				{0.971180894e-3, 80, 300},
				{0.120241800e-2, 3, 500}},
			tol: 1e-11,
		},

		{
			name: "testSpecificVolumePT",
			f:    if97.SpecificVolumePT,
			values: [][]float64{
				{0.394913866e2, 0.0035, 300}, // region 2
				{0.923015898e2, 0.0035, 700},
				{0.542946619e-2, 30, 700}},
			tol: 1e-4,
		},
		{
			name: "testSpecificVolumePT",
			f:    if97.SpecificVolumePT,
			values: [][]float64{
				{1.470853100e-3, 50, 630}, // region 3a
				{1.503831359e-3, 80, 670},
				{2.204728587e-3, 50, 710}, // region 3b
				{1.973692940e-3, 80, 750},
				{1.761696406e-3, 20, 630}, // region 3c
				{1.819560617e-3, 30, 650},
				{2.245587720e-3, 26, 656}, // region 3d
				{2.506897702e-3, 30, 670},
				{2.970225962e-3, 26, 661}, // region 3e
				{3.004627086e-3, 30, 675},
				{5.019029401e-3, 26, 671}, // region 3f
				{4.656470142e-3, 30, 690},
				{2.163198378e-3, 23.6, 649}, // region 3g
				{2.166044161e-3, 24, 650},
				{2.651081407e-3, 23.6, 652}, // region 3h
				{2.967802335e-3, 24, 654},
				{3.273916816e-3, 23.6, 653}, // region 3i
				{3.550329864e-3, 24, 655},
				{4.545001142e-3, 23.5, 655}, // region 3j
				{5.100267704e-3, 24, 660},
				{6.109525997e-3, 23, 660}, // region 3k
				{6.427325645e-3, 24, 670},
				{2.117860851e-3, 22.6, 646}, // region 3l
				{2.062374674e-3, 23, 646},
				{2.533063780e-3, 22.6, 648.6}, // region 3m
				{2.572971781e-3, 22.8, 649.3},
				{2.923432711e-3, 22.6, 649.0}, // region 3n
				{2.913311494e-3, 22.8, 649.7},
				{3.131208996e-3, 22.6, 649.1}, // region 3o
				{3.221160278e-3, 22.8, 649.9},
				{3.715596186e-3, 22.6, 649.4}, // region 3p
				{3.664754790e-3, 22.8, 650.2},
				{1.970999272e-3, 21.1, 640}, // region 3q
				{2.043919161e-3, 21.8, 643},
				{5.251009921e-3, 21.1, 644}, // region 3r
				{5.256844741e-3, 21.8, 648},
				{1.932829079e-3, 19.1, 635}, // region 3s
				{1.985387227e-3, 20, 638},
				{8.483262001e-3, 17, 626}, // region 3t
				{6.227528101e-3, 20, 640}},
			tol: 1e-12,
		},

		{
			name: "testSpecificVolumePT",
			f:    if97.SpecificVolumePT,
			values: [][]float64{
				{0.138455090e1, 0.5, 1500}, // region 5
				{0.230761299e-1, 30, 1500},
				{0.311385219e-1, 30, 2000}},
			tol: 1e-8,
		},

		{
			name: "testSpecificVolumePT",
			f:    if97.SpecificVolumePT,
			values: [][]float64{
				{0.138455090e1, 0.5, 1500}, // region 5
				{0.230761299e-1, 30, 1500},
				{0.311385219e-1, 30, 2000}},
			tol: 1e-8,
		},

		{
			name: "testSpeedOfSoundPT",
			f:    if97.SpeedOfSoundPT,
			values: [][]float64{
				{0.150773921e4, 3, 300}, // region 1
				{0.163469054e4, 80, 300},
				{0.124071337e4, 3, 500},
				{0.427920172e3, 0.0035, 300}, // region 2
				{0.644289068e3, 0.0035, 700},
				{0.480386523e3, 30, 700},
				//{0.498408101e3, 1, 450}, // region 2 meta
				//{0.489363295e3, 1, 440},
				//{0.481941819e3, 1.5, 450},
				{0.917068690e3, 0.5, 1500}, // region 5
				{0.928548002e3, 30, 1500},
				{0.106736948e4, 30, 2000}},
			tol: 1e-5,
		},

		{
			name: "testSpeedOfSoundPT",
			f:    if97.SpeedOfSoundPT,
			values: [][]float64{
				{0.502005554e3, thirdRegion.REGION3.PressureRhoT(500, 650), 650}, // region 3
				{0.383444594e3, thirdRegion.REGION3.PressureRhoT(200, 650), 650},
				{0.760696041e3, thirdRegion.REGION3.PressureRhoT(500, 750), 750}},
			tol: 1e-2,
		},

		{
			name: "testTemperatureHS",
			f:    if97.TemperatureHS,
			values: [][]float64{
				{3.468475498e2, 1800, 5.3}, // region 4
				{4.251373305e2, 2400, 6},
				{5.225579013e2, 2500, 5.5}},
			tol: 1e-5,
		},

		{
			name: "testTemperaturePH",
			f:    if97.TemperaturePH,
			values: [][]float64{
				{0.391798509e3, 3, 500}, // region 1
				{0.378108626e3, 80, 500},
				{0.611041229e3, 80, 1500},
				{0.534433241e3, 0.001, 3000}, // region 2a
				{0.575373370e3, 3, 3000},
				{0.101077577e4, 3, 4000},
				{0.801299102e3, 5, 3500}, // region 2b
				{0.101531583e4, 5, 4000},
				{0.875279054e3, 25, 3500},
				{0.743056411e3, 40, 2700}, // region 2c
				{0.791137067e3, 60, 2700},
				{0.882756860e3, 60, 3200},
				{6.293083892e2, 20, 1700}, // region 3
				{6.905718338e2, 50, 2000},
				{7.336163014e2, 100, 2100},
				{6.418418053e2, 20, 2500},
				{7.351848618e2, 50, 2400},
				{8.420460876e2, 100, 2700}},
			tol: 1e-5,
		},

		{
			name: "testTemperaturePS",
			f:    if97.TemperaturePS,
			values: [][]float64{
				{0.307842258e3, 3, 0.5}, // region 1
				{0.309979785e3, 80, 0.5},
				{0.565899909e3, 80, 3},
				{0.399517097e3, 0.1, 7.5}, // region 2
				{0.514127081e3, 0.1, 8},
				{0.103984917e4, 2.5, 8},
				{0.600484040e3, 8, 6},
				{0.106495556e4, 8, 7.5},
				{0.103801126e4, 90, 6},
				{0.697992849e3, 20, 5.75},
				{0.854011484e3, 80, 5.25},
				{0.949017998e3, 80, 5.75},
				{6.282959869e2, 20, 3.8}, // region 3
				{6.297158726e2, 50, 3.6},
				{7.056880237e2, 100, 4.0},
				{6.401176443e2, 20, 5.0},
				{7.163687517e2, 50, 4.5},
				{8.474332825e2, 100, 5.0}},
			tol: 1e-5,
		},
		{
			name  : "testDielectricConstantPT",
			f:if97.DielectricConstantPT,
			values :[][]float64{
				{0.785907250e2, 5, 298.15},
				//meta {0.112620970e1, 10, 873.15},
				{0.103126058e2, 40, 673.15}},
			tol:1e-4,
		},
	
		{
			name  : "testCompressionFactorPT",
			f:if97.CompressionFactorPT,
			values :[][]float64{
				{0.00073, 0.1, 298.15},
				{0.03636, 5, 298.15},
				{0.24600, 40, 673.15}},
			tol:1e-5,
			},

			{
				name  : "testDynamicViscosityPT",
				f:if97.DynamicViscosityPT,
				values :[][]float64{
					{0.890022551e-3, .1, 298.15},
					{0.339743835e-4, 20, 873.15},
					{0.726093560e-4, 60, 673.15}},
				tol:1e-12,
				},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			for _, x := range tc.values {
				//fmt.Println("_____________________",if97.Region.GetName())
				retval, err := tc.f(x[1], x[2])
				assert.Equal(t, err, nil)

				assert.InDelta(t, x[0], retval, tc.tol)
			}
		})
	}



	cases1f := []struct {
		f      func(p float64)  (float64, error)
		name   string
		values [][]float64
		tol    float64
	}{
		{
			f    :  if97.SaturationPressureT,
			name  : "testSaturationPressureT",
			values :[][]float64{
				{0.353658941e-2, 300},
				{0.263889776e1, 500},
				{0.123443146e2, 600}},
			tol  :  1e-7,
		},
		{
			f    :  if97.SaturationTemperatureP,
			name  : "testSaturationTemperatureP",
			values :[][]float64{
				{0.372755919e3, 0.1},
				{0.453035632e3, 1},
				{0.584149488e3, 10}},
			tol  :  1e-6,
		},
		{
			f    :  if97.SurfaceTensionT,
			name  : "testSurfaceTensionT",
			values :[][]float64{
				{0.0716859625, 300},
				{0.0428914992, 450},
				{0.00837561087, 600},
			},
			tol  :  1e-10,
		},
	}
	for _, tc := range cases1f {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			for _, x := range tc.values {
				//fmt.Println("_____________________",if97.Region.GetName())
				retval, err := tc.f(x[1])
				assert.Equal(t, err, nil)

				assert.InDelta(t, x[0], retval, tc.tol)
			}
		})
	}


	t.Run("testThermalConductivityRhoT", func(t *testing.T) {
		t.Parallel()
		for _, x := range [][]float64{
	        {0.607509806, 0.1, 298.15},
	        //{0.867570353e-1, 10, 873.15},
	        {0.398506911, 40, 673.15}} {
			//fmt.Println("_____________________",if97.Region.GetName())
			dens,err := if97.DensityPT(x[1], x[2])
			assert.Equal(t, err, nil)
			retval, err := if97.ThermalConductivityRhoT(dens,x[2])
			assert.Equal(t, err, nil)

			assert.InDelta(t, x[0], retval, 1e-6)
		}
	})

	t.Run("testRefractiveIndexPTLambda", func(t *testing.T) {
		t.Parallel()
		for _, x := range [][]float64{
			{0.139277824e1, 0.1, 298.15, 0.2265},
			{0.133285819e1, 0.1, 298.15, 0.5893},
			//meta {0.101098988e1, 10, 773.15, 0.2265},
			//meta {0.100949307e1, 10, 773.15, 0.5893},
			{0.119757252e1, 40, 673.15, 0.2265},
			{0.116968699e1, 40, 673.15, 0.5893}} {
			//fmt.Println("_____________________",if97.Region.GetName())

			retval, err := if97.RefractiveIndexPTLambda(x[1],x[2],x[3])
			assert.Equal(t, err, nil)

			assert.InDelta(t, x[0], retval, 1e-6)
		}
	})

}
