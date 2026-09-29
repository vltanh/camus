package score

import (
	"math"
	"testing"

	gr "github.com/jsdoublel/camus/internal/graphs"
)

func TestFractionalScoresEqualScaledCounts(t *testing.T) {
	td := makeTreeDataWithQuartets(t, "(((A,B)a,(C,D)b)e,(E,(F,G)f)c)r;", []quartetCount{
		{nwk: "((A,E),(B,F));", count: 7},
		{nwk: "((A,F),(B,E));", count: 3},
	})
	weights := map[gr.Quartet]float64{}
	for _, q := range td.Quartets(td.Root().Id()) {
		weights[q] = float64(td.NumQuartet(q)) / 8
	}
	weighted := gr.MakeWeightedTreeData(&td.Tree, weights)
	var countScorer MaximizeScorer
	if err := countScorer.Init(td, 2); err != nil {
		t.Fatal(err)
	}
	var weightScorer WeightedMaximizeScorer
	if err := weightScorer.Init(weighted, 2); err != nil {
		t.Fatal(err)
	}
	positive := 0
	for u := range td.Nodes() {
		for w := range td.Nodes() {
			want := float64(countScorer.CalcScore(u, w, td)) / 8
			got := weightScorer.CalcScore(u, w, weighted)
			if got != want {
				t.Fatalf("edge %d->%d weighted score=%g, want=%g", u, w, got, want)
			}
			if got > 0 {
				positive++
				branches := []gr.Branch{{IDs: [2]int{u, w}}}
				percent, err := weightScorer.PercentQuartetSat(branches, weighted)
				if err != nil || math.Abs(percent-100*want/1.25) > 1e-12 {
					t.Fatalf("weighted percentage=%g, error=%v", percent, err)
				}
			}
		}
	}
	if positive == 0 {
		t.Fatal("test has no positive candidate scores")
	}
	if err := weightScorer.Init(weighted, 1, AsSet(true)); err == nil {
		t.Fatal("accepted -asSet with weights")
	}
	if err := weightScorer.Init(td, 1); err == nil {
		t.Fatal("accepted missing weights")
	}
}
