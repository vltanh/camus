// Network/CSV compatibility serializer for the C++ inference engine.
// This process never reads quartet counts or evidence weights.
package main

import (
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/evolbioinfo/gotree/io/newick"
	gr "github.com/jsdoublel/camus/internal/graphs"
	pr "github.com/jsdoublel/camus/internal/prep"
)

type result struct {
	Tree       string
	Signed     bool
	Candidates []struct {
		Score    string
		Branches [][2]int
	}
}

func run() error {
	if len(os.Args) != 3 {
		return fmt.Errorf("usage: camus-format branches.json output.csv")
	}
	data, err := os.ReadFile(os.Args[1])
	if err != nil {
		return err
	}
	var input result
	if err = json.Unmarshal(data, &input); err != nil {
		return err
	}
	text, err := os.ReadFile(input.Tree)
	if err != nil {
		return err
	}
	tre, err := newick.NewParser(strings.NewReader(string(text))).Parse()
	if err != nil {
		return err
	}
	tre.ClearLengths(true, true)
	tre.ClearComments()
	tre.ClearSupports()
	tre.RemoveSingleNodes()
	for i, n := range tre.Nodes() {
		n.SetId(i)
	}
	if err = tre.UpdateTipIndex(); err != nil {
		return err
	}
	td := gr.MakeTreeData(tre, nil)
	newicks := make([]string, len(input.Candidates))
	scores := make([]float64, len(input.Candidates))
	for i, candidate := range input.Candidates {
		if len(candidate.Branches) != i+1 {
			return fmt.Errorf("incorrect branch count for candidate %d", i+1)
		}
		branches := make([]gr.Branch, len(candidate.Branches))
		for j, ids := range candidate.Branches {
			if ids[0] <= 0 || ids[1] <= 0 || ids[0] >= len(td.Nodes()) || ids[1] >= len(td.Nodes()) {
				return fmt.Errorf("branch node IDs out of range")
			}
			branches[j] = gr.Branch{IDs: ids}
		}
		newicks[i] = gr.MakeNetwork(td, branches).Newick()
		if scores[i], err = strconv.ParseFloat(candidate.Score, 64); err != nil {
			return err
		}
	}
	write := pr.WriteDPResultsToCSV
	if input.Signed {
		write = pr.WriteSignedDPResultsToCSV
	}
	file, err := os.Create(os.Args[2])
	if err != nil {
		return err
	}
	defer file.Close()
	if err = write(td, newicks, scores, file); err != nil {
		return err
	}
	if err = write(td, newicks, scores, os.Stdout); err != nil {
		return err
	}
	if !input.Signed {
		return pr.WriteResultsLineplot(scores, strings.TrimSuffix(os.Args[2], ".csv"))
	}
	return nil
}

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}
