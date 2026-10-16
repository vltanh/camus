package score

import (
	"math"
	"strings"
	"testing"

	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
	gr "github.com/jsdoublel/camus/internal/graphs"
)

func exhaustiveEvidence(t testing.TB, nwk string) *gr.TreeData {
	t.Helper()
	tre, err := newick.NewParser(strings.NewReader(nwk)).Parse()
	if err != nil {
		t.Fatal(err)
	}
	if err := tre.UpdateTipIndex(); err != nil {
		t.Fatal(err)
	}
	for i, node := range tre.Nodes() {
		node.SetId(i)
	}
	weights := map[gr.Quartet]float64{}
	n := len(tre.AllTipNames())
	for d := 3; d < n; d++ {
		for c := 2; c < d; c++ {
			for b := 1; b < c; b++ {
				for a := 0; a < b; a++ {
					for i, ids := range [][4]int16{{int16(a), int16(b), int16(c), int16(d)},
						{int16(a), int16(c), int16(b), int16(d)}, {int16(a), int16(d), int16(b), int16(c)}} {
						q := gr.QuartetFromTreeQ(&tree.Quartet{T1: 0, T2: 1, T3: 2, T4: 3}, ids[:])
						// Signed, non-dyadic weights with a broad magnitude range expose
						// changes to the floating-point summation order.
						weights[q] = math.Ldexp(float64((a+3*b+7*c+11*d+i)%19-9)/13, (a+b+c+d+i)%23-11)
					}
				}
			}
		}
	}
	return gr.MakeWeightedTreeData(tre, weights)
}

func TestCachedEdgeCompatibilityAndWeightedScoresAreBitExact(t *testing.T) {
	for _, nwk := range []string{
		"(OUT,(((A,B),(C,D)),(E,(F,G))));",
		"(((((((OUT,A),B),C),D),E),F),G);",
		"((OUT,A),((B,(C,D)),(E,(F,G))));",
	} {
		td := exhaustiveEvidence(t, nwk)
		var weighted WeightedMaximizeScorer
		if err := weighted.Init(td, 2); err != nil {
			t.Fatal(err)
		}
		matched := 0
		for u := range td.Nodes() {
			for w := range td.Nodes() {
				if !ShouldCalcEdge(u, w, td) {
					continue
				}
				v := td.LCA(u, w)
				edge := newEdgeQuartetScorer(u, w, td)
				wSub := getWSubtree(u, w, v, td)
				var expected float64
				for _, q := range td.Quartets(v) {
					want := QuartetScore(q, td.IdToNodes[u], td.IdToNodes[w], td.IdToNodes[v], wSub, td) == gr.Qeq
					if got := edge.matches(q); got != want {
						t.Fatalf("edge %d->%d quartet %x cached=%v reference=%v", u, w, q, got, want)
					}
					if want {
						expected += td.QuartetWeight(q)
						matched++
					}
				}
				if got := weighted.CalcScore(u, w, td); math.Float64bits(got) != math.Float64bits(expected) {
					t.Fatalf("edge %d->%d changed weighted sum: %.17g vs %.17g", u, w, got, expected)
				}
			}
		}
		if matched == 0 {
			t.Fatal("fixture has no matching quartets")
		}
	}
}

var quartetBenchResult int

func BenchmarkEdgeQuartetCompatibility(b *testing.B) {
	td := exhaustiveEvidence(b, "(OUT,(((A,B),(C,D)),(E,(F,G))));")
	for _, cached := range []bool{false, true} {
		name := "reference"
		if cached {
			name = "cached"
		}
		b.Run(name, func(b *testing.B) {
			matches := 0
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				for u := range td.Nodes() {
					for w := range td.Nodes() {
						if !ShouldCalcEdge(u, w, td) {
							continue
						}
						v := td.LCA(u, w)
						if cached {
							edge := newEdgeQuartetScorer(u, w, td)
							for _, q := range td.Quartets(v) {
								if edge.matches(q) {
									matches++
								}
							}
						} else {
							wSub := getWSubtree(u, w, v, td)
							for _, q := range td.Quartets(v) {
								if QuartetScore(q, td.IdToNodes[u], td.IdToNodes[w], td.IdToNodes[v], wSub, td) == gr.Qeq {
									matches++
								}
							}
						}
					}
				}
			}
			quartetBenchResult = matches
		})
	}
}
