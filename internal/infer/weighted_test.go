package infer

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"testing"

	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
	sc "github.com/jsdoublel/camus/internal/score"
)

// Fractional weights below one must influence the final network, including when
// the preferred conflicting quartet changes. Scaling evidence must preserve it.
func TestWeightedInferenceMatchesScaledCounts(t *testing.T) {
	directory := t.TempDir()
	treePath := filepath.Join(directory, "constraint.nwk")
	if err := os.WriteFile(treePath, []byte("((A,B),(C,D));\n"), 0600); err != nil {
		t.Fatal(err)
	}
	run := func(first, second float64, weighted bool) *DPResults {
		t.Helper()
		field, format := "count", pr.QuartetCounts
		if weighted {
			field, format = "weight", pr.QuartetWeights
		}
		path := filepath.Join(directory, field+".csv")
		rows := fmt.Sprintf("a,b,c,d,%s\nA,B,C,D,100\nA,C,B,D,%.17g\nA,D,B,C,%.17g\n", field, first, second)
		if err := os.WriteFile(path, []byte(rows), 0600); err != nil {
			t.Fatal(err)
		}
		tre, genes, err := pr.ReadInputFiles(treePath, path, format)
		if err != nil {
			t.Fatal(err)
		}
		opts := BuildTestInferOpts(t, 0, 0.5, &sc.MaximizeScorer{}, 0)
		opts.NProcs = 1
		if weighted {
			opts.QuartetOpts.WeightsFile = path
		} else {
			opts.QuartetOpts.CountsFile = path
		}
		results, err := Infer(tre, genes.Trees, opts)
		if err != nil {
			t.Fatal(err)
		}
		return results
	}
	var winners []string
	for _, pair := range [][2]float64{{6, 1}, {1, 6}} {
		counts := run(pair[0], pair[1], false)
		if len(counts.Branches) != 1 {
			t.Fatalf("expected one reticulation for conflicting four-taxon evidence, got %d", len(counts.Branches))
		}
		expected := gr.MakeNetwork(counts.Tree, counts.Branches[0]).Newick()
		winners = append(winners, expected)
		for _, factor := range []float64{1, 0.125, 1e-12} {
			weights := run(pair[0]*factor, pair[1]*factor, true)
			if len(weights.Branches) != len(counts.Branches) {
				t.Fatalf("fractional weights changed result count: %d vs %d", len(weights.Branches), len(counts.Branches))
			}
			actual := gr.MakeNetwork(weights.Tree, weights.Branches[0]).Newick()
			if actual != expected {
				t.Fatalf("scaled weights chose %s, counts chose %s", actual, expected)
			}
			if math.Abs(weights.QSatScore[0]-100*6.0/7) > 1e-10 {
				t.Fatalf("incorrect weighted satisfaction=%g", weights.QSatScore[0])
			}
		}
	}
	if winners[0] == winners[1] {
		t.Fatal("changing the preferred quartet did not change the winning network")
	}
	zero := run(0, 0, true)
	if len(zero.Branches) != 0 || zero.Tree.TotalQuartetWeight() != 0 {
		t.Fatal("zero evidence should return the constraint tree without reticulations")
	}
}
