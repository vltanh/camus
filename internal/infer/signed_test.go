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

func TestSignedSearchUsesNegativeEvidenceAndRawGain(t *testing.T) {
	dir := t.TempDir()
	treePath := filepath.Join(dir, "tree.nwk")
	if err := os.WriteFile(treePath, []byte("((A,B),(C,D));\n"), 0600); err != nil {
		t.Fatal(err)
	}
	for _, tc := range []struct {
		first, second float64
		edges         int
		gain          float64
	}{
		{0.25, -0.9, 1, 0.25}, {-0.25, 0.4, 1, 0.4}, {-0.25, -0.4, 0, 0}, {0, 0, 0, 0},
	} {
		for _, backbone := range []float64{-1000, 1000} {
			path := filepath.Join(dir, "weights.csv")
			rows := fmt.Sprintf("a,b,c,d,weight\nA,B,C,D,%.17g\nA,C,B,D,%.17g\nA,D,B,C,%.17g\n", backbone, tc.first, tc.second)
			if err := os.WriteFile(path, []byte(rows), 0600); err != nil {
				t.Fatal(err)
			}
			tre, genes, err := pr.ReadInputFiles(treePath, path, pr.QuartetWeights)
			if err != nil {
				t.Fatal(err)
			}
			opts := BuildTestInferOpts(t, 0, 0.5, &sc.MaximizeScorer{}, 0)
			opts.NProcs = 1
			opts.QuartetOpts.WeightsFile = path
			opts.QuartetOpts.SignedWeights = true
			opts.QuartetOpts.MaxReticulations = -1
			result, err := Infer(tre, genes.Trees, opts)
			if err != nil {
				t.Fatal(err)
			}
			bestGain, bestEdges := 0., 0
			for i, gain := range result.QSatScore {
				if len(result.Branches[i]) != i+1 {
					t.Fatalf("phantom cardinality %d: %v", i+1, result.Branches[i])
				}
				if gain > bestGain {
					bestGain, bestEdges = gain, len(result.Branches[i])
				}
			}
			if bestEdges != tc.edges || math.Abs(bestGain-tc.gain) > 1e-12 {
				t.Fatalf("weights %g,%g: edges=%d gain=%g, want %d,%g", tc.first, tc.second, bestEdges, bestGain, tc.edges, tc.gain)
			}
		}
	}
}

// One long cycle (gain 10) blocks three mutually compatible short cycles
// (gain 4 each). Exact k=2 has gain 8: stopping there misses k=3, gain 12.
type plateauScorer struct{}

func (plateauScorer) Init(*gr.TreeData, int, ...sc.ScoreOptions) error { return nil }
func (plateauScorer) CalcScore(u, w int, td *gr.TreeData) float64 {
	from, to := td.IdToNodes[u].Name(), td.IdToNodes[w].Name()
	if from == "Z" && to == "L" {
		return 10
	}
	if (from == "A" && to == "C") || (from == "E" && to == "G") || (from == "I" && to == "K") {
		return 4
	}
	return -100
}
func (s plateauScorer) PercentQuartetSat(branches []gr.Branch, td *gr.TreeData) (float64, error) {
	var gain float64
	for _, branch := range branches {
		gain += s.CalcScore(branch.IDs[0], branch.IDs[1], td)
	}
	return gain, nil
}

func TestSignedCapContinuesAfterCardinalityPlateau(t *testing.T) {
	dir := t.TempDir()
	treePath := filepath.Join(dir, "tree.nwk")
	countPath := filepath.Join(dir, "empty.csv")
	if err := os.WriteFile(treePath, []byte("(Z,(A,(B,(C,(D,(E,(F,(G,(H,(I,(J,(K,L))))))))))));\n"), 0600); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(countPath, []byte("a,b,c,d,count\n"), 0600); err != nil {
		t.Fatal(err)
	}
	tre, genes, err := pr.ReadInputFiles(treePath, countPath, pr.QuartetCounts)
	if err != nil {
		t.Fatal(err)
	}
	opts, err := pr.SetQuartetFilterOptions(0, .5)
	if err != nil {
		t.Fatal(err)
	}
	opts.CountsFile = countPath
	td, err := pr.Preprocess(tre, genes.Trees, 1, opts, 0)
	if err != nil {
		t.Fatal(err)
	}
	dp, err := newDP[float64](plateauScorer{}, td, 1)
	if err != nil {
		t.Fatal(err)
	}
	dp.Capped = true
	dp.MaxEdges = dp.NumNodes // Uncapped: only the constraint tree's structural bound.
	result := dp.RunDP()
	if len(result.QSatScore) != 3 {
		t.Fatalf("lost cardinality states: %v", result.QSatScore)
	}
	for i, want := range []float64{10, 8, 12} {
		if result.QSatScore[i] != want {
			t.Fatalf("k=%d: gain %g, want %g", i+1, result.QSatScore[i], want)
		}
	}
}
