package region

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestRegion(t *testing.T) {

	cases := []struct {
		name   string
		f      func(arg float64) float64
		values [][]float64
		tol    float64
	}{
		{
		name  :"testPressureB23",
		f  :    pressureB23,
		values: [][]float64{{0.1652916425e2, 0.623150000e3}},
		tol:    1e-8,
		},
		{
		name  :"testSaturationPressure3",
		f  :    saturationPressure3,
		values: [][]float64{
			{1.687755057e1, 3.8},
            {2.164451789e1, 4.2},
            {1.668968482e1, 5.2},},
			tol:    1e-8,
		},	
		{
		name  :"testSpecificEnthalpy",
		f  :    specificEnthalpy1,
		values: [][]float64{
			{3.085509647e2, 1.0},
			{7.006304472e2, 2.0},
			{1.198359754e3, 3.0}},
		tol:    1e-6,
		},
		{
		name  :"testSpecificEnthalpy2ab",
		f  :    specificEnthalpy2ab,
		values: [][]float64{
			{2.723729985e3, 7.0},
			{2.599047210e3, 8.0},
			{2.511861477e3, 9.0}},
		tol:    1e-6,
		},
		{
		name  :"testSpecificEnthalpy2ab",
		f  :    specificEnthalpy2ab,
		values: [][]float64{
			{2.723729985e3, 7.0},
			{2.599047210e3, 8.0},
			{2.511861477e3, 9.0}},
			tol:    1e-6,
		},
		{
		name  :"testSpecificEnthalpy2c3b",
		f  :    specificEnthalpy2c3b,
		values: [][]float64{ 
			{2.687693850e3, 5.5},
			{2.451623609e3, 5.0},
			{2.144360448e3, 4.5}},
		tol:    1e-6,
		},
		{
		name  :"testSpecificEnthalpy3a",
		f  :    specificEnthalpy3a,
		values: [][]float64{
			{1.685025565e3, 3.8},
			{1.816891476e3, 4.0},
			{1.949352563e3, 4.2}},
		tol:    1e-6,
		},
		{
		name  :"testSpecificEnthalpyB13",
		f  :    specificEnthalpyB13,
		values: [][]float64{
			{1.632525047e3, 3.7},
			{1.593027214e3, 3.6},
			{1.566104611e3, 3.5}},
		tol:    1e-6,
		},
		{
		name  :"testTemperatureB23P",
		f  :    temperatureB23P,
		values: [][]float64{
			{0.623150000e3, 0.1652916425e2}},
		tol:    1e-7,
		},

			
		
	}
		for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			for _, x := range tc.values {
				assert.InDelta(t, x[0], tc.f(x[1]), tc.tol)
			}
		})
	}




	var IF97test = IF97Region{}
	selectCases:= []struct {
		name   string
		f      func(pressure float64, temperature float64) (int,error)
		args [][]float64
		expected  []int
		tol    float64
	}{
		{
			name:"testRegionHS",
		f  :    IF97test.GetRegionHS,
		args: [][]float64{
			{1200, 3},
			{3000, 7},
			{1900, 4},
			{1000, 3},
		},
		expected:  []int{1, 2, 3, 4},

		},
		{
			name:"testRegionPH",
		f  :    IF97test.GetRegionPH,
		args: [][]float64{
			{10, 1000},
            {10, 3000},
            {25, 2000},
            {10, 2000},
		},
		expected:  []int{ 1, 2, 3, 4},

		},
	}
	for _, tc := range selectCases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			for i, x := range tc.args {
				region,err := tc.f(x[0],x[1]);
				if err != nil{
					fmt.Printf("FAIL TEST %s/n",tc.name)
				}
				assert.Equal(t, tc.expected[i], region)
			}
		})



	}
	t.Run("testTemperatureB23HS", func(t *testing.T) {
		t.Parallel()
		for _, x := range [][]float64{
			{7.135259364e2, 2600, 5.1},
			{7.685345532e2, 2700, 5.15},
			{8.176202120e2, 2800, 5.2},
		}{ 
		assert.InDelta(t,x[0], temperatureB23HS(x[1], x[2]), 1e-7);
	}})


}