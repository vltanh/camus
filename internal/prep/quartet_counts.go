package prep

import (
	"encoding/csv"
	"fmt"
	"github.com/evolbioinfo/gotree/tree"
	gr "github.com/jsdoublel/camus/internal/graphs"
	"io"
	"os"
	"strconv"
)

// Each CSV row a,b,c,d,count represents the unrooted topology a,b|c,d.
// Equivalent rows are canonicalized and their multiplicities summed.
func readQuartetCounts(path string, tre *tree.Tree) (map[gr.Quartet]uint32, error) {
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
		if _, ok := fields[h]; ok {
			return nil, fmt.Errorf("duplicate CSV column %q", h)
		}
		fields[h] = i
	}
	for _, key := range []string{"a", "b", "c", "d", "count"} {
		if _, ok := fields[key]; !ok {
			return nil, fmt.Errorf("missing CSV column %s", key)
		}
	}
	counts := map[gr.Quartet]uint32{}
	for line := 2; ; line++ {
		row, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
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
		c, err := strconv.ParseUint(row[fields["count"]], 10, 32)
		if err != nil {
			return nil, fmt.Errorf("line %d: invalid count: %w", line, err)
		}
		q := gr.QuartetFromTreeQ(&tree.Quartet{T1: 0, T2: 1, T3: 2, T4: 3}, ids)
		total := uint64(counts[q]) + c
		if total > uint64(^uint32(0)) {
			return nil, fmt.Errorf("line %d: accumulated count exceeds uint32", line)
		}
		if total > 0 {
			counts[q] = uint32(total)
		}

	}
	return counts, nil
}
