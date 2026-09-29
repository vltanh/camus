package prep

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
)

func TestQuartetWeightsInput(t *testing.T) {
	tre, err := newick.NewParser(strings.NewReader("((a,b),(c,d));")).Parse()
	if err != nil {
		t.Fatal(err)
	}
	if err := tre.UpdateTipIndex(); err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		name, rows string
		want       float64
		fail       bool
	}{
		{name: "fractional", rows: "a,b,c,d,0.125\na,c,b,d,0.75\n", want: 0.875},
		{name: "equivalent", rows: "a,b,c,d,0.125\nd,c,b,a,0.25\n", want: 0.375},
		{name: "tiny", rows: "a,b,c,d,1e-100\n", want: 1e-100},
		{name: "zero", rows: "a,b,c,d,0\n"},
		{name: "empty"},
		{name: "unknown", rows: "a,b,c,x,0.5\n", fail: true},
		{name: "repeated taxon", rows: "a,b,c,c,0.5\n", fail: true},
		{name: "negative", rows: "a,b,c,d,-0.5\n", fail: true},
		{name: "nan", rows: "a,b,c,d,NaN\n", fail: true},
		{name: "infinite", rows: "a,b,c,d,+Inf\n", fail: true},
		{name: "overflow", rows: "a,b,c,d,1e309\n", fail: true},
		{name: "duplicate overflow", rows: "a,b,c,d,1e308\nd,c,b,a,1e308\n", fail: true},
		{name: "total overflow", rows: "a,b,c,d,1e308\na,c,b,d,1e308\n", fail: true},
		{name: "short row", rows: "a,b,c,d\n", fail: true},
		{name: "not a number", rows: "a,b,c,d,x\n", fail: true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "weights.csv")
			if err := os.WriteFile(path, []byte("a,b,c,d,weight\n"+tc.rows), 0600); err != nil {
				t.Fatal(err)
			}
			weights, err := readQuartetWeights(path, tre)
			if (err != nil) != tc.fail {
				t.Fatalf("error=%v, expected failure=%v", err, tc.fail)
			}
			if tc.fail {
				return
			}
			var sum float64
			for _, weight := range weights {
				if weight <= 0 || math.IsNaN(weight) || math.IsInf(weight, 0) {
					t.Fatalf("invalid retained weight: %g", weight)
				}
				sum += weight
			}
			if sum != tc.want {
				t.Fatalf("weight sum %.17g, want %.17g", sum, tc.want)
			}
			if tc.name == "equivalent" && len(weights) != 1 {
				t.Fatalf("equivalent quartets not merged: %v", weights)
			}
		})
	}
	for _, header := range []string{"a,b,c,d,count", "a,b,c,d,d", "a,b,c,weight"} {
		path := filepath.Join(t.TempDir(), "weights.csv")
		if err := os.WriteFile(path, []byte(header+"\n"), 0600); err != nil {
			t.Fatal(err)
		}
		if _, err := readQuartetWeights(path, tre); err == nil {
			t.Fatalf("accepted invalid header %q", header)
		}
	}
}

func TestWeightedPreprocessingPreservesEvidence(t *testing.T) {
	path := filepath.Join(t.TempDir(), "weights.csv")
	if err := os.WriteFile(path, []byte("a,b,c,d,weight\na,b,c,d,10\na,c,b,d,0.125\nd,b,c,a,0.25\na,d,b,c,0\n"), 0600); err != nil {
		t.Fatal(err)
	}
	tre, err := newick.NewParser(strings.NewReader("((a,b),(c,d));")).Parse()
	if err != nil {
		t.Fatal(err)
	}
	opts, err := SetQuartetFilterOptions(0, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	opts.WeightsFile = path
	td, err := Preprocess(tre, nil, 1, opts, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !td.HasQuartetWeights() || td.TotalNumUniqueQuartets() != 1 || td.TotalQuartetWeight() != 0.375 {
		t.Fatalf("incorrect evidence after removing backbone quartets: unique=%d weight=%g", td.TotalNumUniqueQuartets(), td.TotalQuartetWeight())
	}
	for _, invalid := range []QuartetFilterOptions{
		{WeightsFile: path, CountsFile: path},
		{WeightsFile: path, mode: 2},
	} {
		if _, err := Preprocess(tre.Clone(), nil, 1, invalid, 0); err == nil {
			t.Fatal("accepted incompatible weight preprocessing")
		}
	}
	if _, err := Preprocess(tre.Clone(), nil, 1, opts, 0.5); err == nil {
		t.Fatal("accepted support collapse with weighted input")
	}
}
