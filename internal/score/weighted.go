package score

import (
	"fmt"
	"log"

	gr "github.com/jsdoublel/camus/internal/graphs"
	"golang.org/x/sync/errgroup"
)

// WeightedMaximizeScorer uses the same quartet compatibility rule as max scoring,
// with fixed floating-point evidence weights in place of occurrence counts.
type WeightedMaximizeScorer struct {
	totals [][]float64
	Signed bool
}

func (s *WeightedMaximizeScorer) Init(td *gr.TreeData, nprocs int, opts ...ScoreOptions) error {
	if !td.HasQuartetWeights() {
		return fmt.Errorf("weighted max scoring requires quartet weights")
	}
	var options scorerOpts
	for _, opt := range opts {
		if err := opt(&options); err != nil {
			return err
		}
	}
	if options.asSet {
		return fmt.Errorf("quartet weights cannot be scored as a set")
	}
	log.Println("calculating weighted edge scores")
	n := len(td.Nodes())
	s.totals = make([][]float64, n)
	var tasks errgroup.Group
	tasks.SetLimit(nprocs)
	for u := range n {
		s.totals[u] = make([]float64, n)
		tasks.Go(func() error {
			for w := range n {
				if !ShouldCalcEdge(u, w, td) {
					continue
				}
				v := td.LCA(u, w)
				wSub := getWSubtree(u, w, v, td)
				for _, q := range td.Quartets(v) {
					if QuartetScore(q, td.IdToNodes[u], td.IdToNodes[w], td.IdToNodes[v], wSub, td) == gr.Qeq {
						s.totals[u][w] += td.QuartetWeight(q)
					}
				}
			}
			return nil
		})
	}
	return tasks.Wait()
}

func (s WeightedMaximizeScorer) CalcScore(u, w int, td *gr.TreeData) float64 {
	return s.totals[u][w]
}

// Percentage is relative to the total input weight after removing quartets
// already displayed by the constraint tree, matching count mode's denominator.
func (s WeightedMaximizeScorer) PercentQuartetSat(branches []gr.Branch, td *gr.TreeData) (float64, error) {
	if s.totals == nil {
		return 0, ErrQuartetsNotInit
	}
	var matched float64
	for _, branch := range branches {
		u, w := branch.IDs[0], branch.IDs[1]
		if u < 0 || w < 0 || u >= len(s.totals) || w >= len(s.totals) {
			return 0, fmt.Errorf("node ids [%d, %d] out of range %d", u, w, len(s.totals))
		}
		matched += s.totals[u][w]
	}
	if s.Signed {
		return matched, nil // Raw gain; a signed-weight sum is not a percentage denominator.
	}
	total := td.TotalQuartetWeight()
	if total == 0 {
		return 0, nil
	}
	return 100 * (matched / total), nil
}
