package prep

import (
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"

	"github.com/evolbioinfo/gotree/tree"
	gr "github.com/jsdoublel/camus/internal/graphs"
)

// Each row gives a fixed, nonnegative weight for topology a,b|c,d.
// Equivalent topology rows are canonicalized and summed without integer rounding.
func readQuartetWeights(path string, tre *tree.Tree) (map[gr.Quartet]float64, error) {
	return readEvidenceWeights(path, tre, false)
}

func readEvidenceWeights(path string, tre *tree.Tree, signed bool) (map[gr.Quartet]float64, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	r := csv.NewReader(f)
	header, err := r.Read()
	if err != nil {
		return nil, err
	}
	fields := map[string]int{}
	for i, h := range header {
		if _, exists := fields[h]; exists {
			return nil, fmt.Errorf("duplicate CSV column %q", h)
		}
		fields[h] = i
	}
	for _, key := range []string{"a", "b", "c", "d", "weight"} {
		if _, exists := fields[key]; !exists {
			return nil, fmt.Errorf("missing CSV column %s", key)
		}
	}
	weights := map[gr.Quartet]float64{}
	for line := 2; ; line++ {
		row, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("line %d: %w", line, err)
		}
		ids := make([]int16, 4)
		seen := map[int]bool{}
		for i, key := range []string{"a", "b", "c", "d"} {
			name := row[fields[key]]
			id, err := tre.TipIndex(name)
			if err != nil {
				return nil, fmt.Errorf("line %d: unknown taxon %q", line, name)
			}
			if id < 0 || id >= 1<<15 || seen[id] {
				return nil, fmt.Errorf("line %d: quartet requires four distinct taxa within CAMUS capacity", line)
			}
			seen[id] = true
			ids[i] = int16(id)
		}
		weight, err := strconv.ParseFloat(strings.TrimSpace(row[fields["weight"]]), 64)
		if err != nil || math.IsNaN(weight) || math.IsInf(weight, 0) || (!signed && weight < 0) {
			return nil, fmt.Errorf("line %d: weight must be a finite nonnegative number", line)
		}
		if weight == 0 {
			continue
		}
		q := gr.QuartetFromTreeQ(&tree.Quartet{T1: 0, T2: 1, T3: 2, T4: 3}, ids)
		total := weights[q] + weight
		if math.IsInf(total, 0) {
			return nil, fmt.Errorf("line %d: accumulated quartet weight is not finite", line)
		}
		weights[q] = total
		if total == 0 {
			delete(weights, q)
		}
	}
	keys := make([]gr.Quartet, 0, len(weights))
	for q := range weights {
		keys = append(keys, q)
	}
	slices.Sort(keys)
	var total float64
	for _, q := range keys {
		total += weights[q]
	}
	if math.IsInf(total, 0) {
		return nil, fmt.Errorf("total quartet weight is not finite")
	}
	return weights, nil
}
