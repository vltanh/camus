package graphs

import (
	"slices"
	"strings"
	"testing"

	"github.com/bits-and-blooms/bitset"
	"github.com/evolbioinfo/gotree/io/newick"
	"github.com/evolbioinfo/gotree/tree"
)

func TestQuartetIndexMatchesIndependentFullScanAndSharesIdenticalLists(t *testing.T) {
	for _, nwk := range []string{"(OUT,(((A,B),(C,D)),(E,(F,G))));", "(OUT,((((((A,B),C),D),E),F),G));"} {
		tre, err := newick.NewParser(strings.NewReader(nwk)).Parse()
		if err != nil {
			t.Fatal(err)
		}
		if err := tre.UpdateTipIndex(); err != nil {
			t.Fatal(err)
		}
		for i, n := range tre.Nodes() {
			n.SetId(i)
		}
		observed, err := QuartetsFromTree(tre.Clone(), tre)
		if err != nil {
			t.Fatal(err)
		}
		counts := map[Quartet]uint32{}
		for q := range observed {
			for _, alternate := range q.AllQuartets() {
				counts[alternate] = 1
			}
		}
		td := MakeTreeData(tre, counts)
		for _, node := range tre.Nodes() {
			want := []Quartet{}
			for q := range counts {
				found := 0
				for i := range 4 {
					if td.leafsets[node.Id()].Test(uint(q.Taxon(i))) {
						found++
					}
				}
				if found >= 3 {
					want = append(want, q)
				}
			}
			slices.Sort(want)
			got := td.Quartets(node.Id())
			if !slices.Equal(got, want) {
				t.Fatalf("node %d index differs from full scan", node.Id())
			}
			if cap(got) != len(got) {
				t.Fatalf("node %d retains unused quartet capacity", node.Id())
			}
		}
		root := td.Quartets(tre.Root().Id())
		shared := false
		for _, child := range td.Children[tre.Root().Id()] {
			if !child.Tip() {
				childSet := td.Quartets(child.Id())
				if len(childSet) == len(root) && &childSet[0] == &root[0] {
					shared = true
				}
			}
		}
		if !shared {
			t.Fatal("identical root/outgroup-side lists were copied")
		}
	}
}

func TestMakeTreeData(t *testing.T) {
	testCases := []struct {
		name        string
		tre         string
		quartets    []string
		lca         map[string][][]string
		leafset     map[string][]string
		quartetSets map[string][]string
	}{
		{
			name:     "basic",
			tre:      "((((A,B)a,C)b,D)c,F)r;",
			quartets: []string{"((A,C),(B,D));"},
			lca: map[string][][]string{
				"a": {
					{"A", "B"},
				},
				"b": {
					{"A", "C"},
					{"B", "C"},
				},
				"c": {
					{"A", "D"},
					{"B", "D"},
					{"C", "D"},
				},
				"r": {
					{"A", "F"},
					{"B", "F"},
					{"C", "F"},
					{"D", "F"},
				},
			},
			leafset: map[string][]string{
				"a": {"A", "B"},
				"b": {"A", "B", "C"},
				"c": {"A", "B", "C", "D"},
				"r": {"A", "B", "C", "D", "F"},
			},
			quartetSets: map[string][]string{
				"a": {},
				"b": {},
				"c": {"((A,C),(B,D));"},
				"r": {"((A,C),(B,D));"},
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.tre)).Parse()
			if err != nil {
				t.Error("invalid newick tree; test is written wrong")
			}
			q := make([]*tree.Tree, 0)
			for _, s := range test.quartets {
				tmp, err := newick.NewParser(strings.NewReader(s)).Parse()
				if err != nil {
					t.Error("invalid newick tree; test is written wrong")
				}
				q = append(q, tmp)
			}
			err = tre.UpdateTipIndex()
			if err != nil {
				t.Error(err)
			}
			qc := makeQCounts(t, q, tre)
			treeData := MakeTreeData(tre, qc)
			lca := treeData.lca
			leafset := treeData.leafsets
			quartetSets := treeData.quartetSet
			nLeaves := len(lca)
			for i := range nLeaves {
				for j := range nLeaves {
					if lca[i][j] != lca[j][i] {
						t.Error("lca structure problem")
					}
				}
			}
			assertLCAEqual(t, lca, test.lca, tre)
			assertLeafsetEqual(t, leafset, test.leafset, tre)
			assertQuartetSetsEqual(t, quartetSets, test.quartetSets, tre)
		})
	}
}

func TestCountLeavesBelow(t *testing.T) {
	testCases := []struct {
		name     string
		newick   string
		expected map[string]uint64
	}{
		{
			name:   "chain",
			newick: "((((A,B)a,C)b,D)c,E)r;",
			expected: map[string]uint64{
				"A": 1,
				"B": 1,
				"C": 1,
				"D": 1,
				"E": 1,
				"a": 2,
				"b": 3,
				"c": 4,
				"r": 5,
			},
		},
		{
			name:   "mixed",
			newick: "((A,(B,C)b)a,(D,E)c)r;",
			expected: map[string]uint64{
				"A": 1,
				"B": 1,
				"C": 1,
				"D": 1,
				"E": 1,
				"a": 3,
				"b": 2,
				"c": 2,
				"r": 5,
			},
		},
	}
	for _, test := range testCases {
		t.Run(test.name, func(t *testing.T) {
			tre, err := newick.NewParser(strings.NewReader(test.newick)).Parse()
			if err != nil {
				t.Fatalf("invalid newick tree: %v", err)
			}
			if err := tre.UpdateTipIndex(); err != nil {
				t.Fatalf("failed to update tip index: %v", err)
			}
			children := children(tre)
			counts := countLeavesBelow(tre, children)
			assertLeavesBelow(t, tre, counts, test.expected)
		})
	}
}

func assertLeavesBelow(t *testing.T, tre *tree.Tree, counts []uint64, expected map[string]uint64) {
	t.Helper()
	for label, want := range expected {
		node := getNode(t, label, tre)
		if got := counts[node.Id()]; got != want {
			t.Fatalf("leaves below %s = %d, want %d", label, got, want)
		}
	}
}

func assertLCAEqual(t *testing.T, lca [][]int, expected map[string][][]string, tre *tree.Tree) {
	t.Helper()
	for label, pairs := range expected {
		for _, pair := range pairs {
			if len(pair) != 2 {
				t.Fatalf("expected LCA pair to have size 2, got %d", len(pair))
			}
			node1 := getNode(t, pair[0], tre)
			node2 := getNode(t, pair[1], tre)
			lcaNode := getNode(t, label, tre)
			if lca[node1.Id()][node2.Id()] != lcaNode.Id() {
				t.Fatalf("lca(%s,%s)=%d, want %d", node1.Name(), node2.Name(), lca[node1.Id()][node2.Id()], lcaNode.Id())
			}
		}
	}
}

func assertLeafsetEqual(t *testing.T, leafset []*bitset.BitSet, expected map[string][]string, tre *tree.Tree) {
	t.Helper()
	for label, leaves := range expected {
		node := getNode(t, label, tre)
		leafsetList := leafset[node.Id()]
		for _, leaf := range leaves {
			id, err := tre.TipIndex(leaf)
			if err != nil {
				t.Fatalf("failed to find tip %q: %v", leaf, err)
			}
			if !leafsetList.Test(uint(id)) {
				t.Fatalf("leafset for %s missing tip %s", label, leaf)
			}
		}
	}
}

func assertQuartetSetsEqual(t *testing.T, got [][]Quartet, expected map[string][]string, tre *tree.Tree) {
	t.Helper()
	for label, quartets := range expected {
		node := getNode(t, label, tre)
		for _, quartetString := range quartets {
			qTree, err := newick.NewParser(strings.NewReader(quartetString)).Parse()
			if err != nil {
				t.Fatalf("cannot parse %s as newick tree: %v", quartetString, err)
			}
			q1, err := NewQuartet(qTree, tre)
			if err != nil {
				t.Fatalf("invalid quartet %s: %v", quartetString, err)
			}
			found := false
			for _, q2 := range got[node.Id()] {
				if q1.Compare(q2) == Qeq {
					found = true
					break
				}
			}
			if !found {
				t.Fatalf("expected quartet %s not found for node %s", quartetString, label)
			}
		}
	}
}

func getNode(t *testing.T, label string, tre *tree.Tree) *tree.Node {
	t.Helper()
	nodeList, err := tre.SelectNodes(label)
	if err != nil {
		t.Fatalf("failed to select node %q: %v", label, err)
	}
	if len(nodeList) != 1 {
		t.Fatalf("expected exactly one node labeled %q, got %d", label, len(nodeList))
	}
	return nodeList[0]
}

func makeQCounts(t *testing.T, qList []*tree.Tree, constTree *tree.Tree) map[Quartet]uint32 {
	t.Helper()
	result := make(map[Quartet]uint32)
	for _, qt := range qList {
		q, err := NewQuartet(qt, constTree)
		if err != nil {
			t.Fatalf("invalid quartet in test data: %v", err)
		}
		result[q] += 1
	}
	return result
}
