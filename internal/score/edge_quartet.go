package score

import gr "github.com/jsdoublel/camus/internal/graphs"

// Cache the topology queries shared by every quartet for a candidate edge.
// The compatibility rule and floating-point summation order stay unchanged.
type edgeTip struct {
	fromU, fromW     int
	inW, inWSub, inU bool
}

type edgeQuartetScorer struct {
	tips    []edgeTip
	depths  []int
	nLeaves int
}

func newEdgeQuartetScorer(u, w int, td *gr.TreeData) edgeQuartetScorer {
	v := td.LCA(u, w)
	wSub := getWSubtree(u, w, v, td).Id()
	s := edgeQuartetScorer{tips: make([]edgeTip, td.NLeaves), depths: td.Depths, nLeaves: td.NLeaves}
	for t := range s.tips {
		tip := uint16(t)
		entry := &s.tips[t]
		entry.inW = td.InLeafset(uint16(w), tip)
		entry.inWSub = td.InLeafset(uint16(wSub), tip)
		entry.inU = td.InLeafset(uint16(u), tip)
		if td.InLeafset(uint16(v), tip) {
			id := td.TipToNodeID(tip)
			entry.fromU = td.LCA(u, id)
			entry.fromW = td.LCA(w, id)
		}
	}
	return s
}

func (s edgeQuartetScorer) matches(q gr.Quartet) bool {
	taxa := [4]uint16{q.Taxon(0), q.Taxon(1), q.Taxon(2), q.Taxon(3)}
	bottom := -1
	for i, tip := range taxa {
		if s.tips[tip].inW {
			if bottom != -1 {
				return false
			}
			bottom = i
		}
	}
	if bottom == -1 {
		return false
	}
	bottomInU := s.tips[taxa[bottom]].inU
	var nodes [4]int
	for i, tip := range taxa {
		entry := s.tips[tip]
		nodes[i] = entry.fromU
		if entry.inWSub || bottomInU {
			nodes[i] = entry.fromW
		}
	}
	if dups(nodes) {
		return false
	}
	minW, maxU := s.nLeaves, -1
	best := -1
	taxaInU := false
	for i, tip := range taxa {
		depth := s.depths[nodes[i]]
		inWSub := s.tips[tip].inWSub
		if !taxaInU && inWSub && depth < minW {
			minW, best = depth, i
		} else if !inWSub && depth > maxU {
			taxaInU = true
			maxU, best = depth, i
		}
	}
	return best >= 0 && best != bottom && ((q.Topology()>>best)&1) == ((q.Topology()>>bottom)&1)
}
