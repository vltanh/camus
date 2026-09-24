package prep

import (
	"github.com/evolbioinfo/gotree/io/newick"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestQuartetCountsInput(t *testing.T) {
	tre, err := newick.NewParser(strings.NewReader("((a,b),(c,d));")).Parse()
	if err != nil {
		t.Fatal(err)
	}
	if err = tre.UpdateTipIndex(); err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		name, rows string
		fail       bool
	}{
		{"valid", "a,b,c,d,7\na,c,b,d,2\na,d,b,c,1\n", false},
		{"zeros", "a,b,c,d,0\n", false},
		{"equivalent", "a,b,c,d,3\nd,c,b,a,4\n", false},
		{"unknown", "a,b,c,x,7\n", true},
		{"repeated taxon", "a,b,c,c,7\n", true},
		{"negative", "a,b,c,d,-1\n", true},
		{"overflow", "a,b,c,d,4294967296\n", true},
		{"sum overflow", "a,b,c,d,4294967295\nd,c,b,a,1\n", true},
		{"short row", "a,b,c,d\n", true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "counts.csv")
			if err := os.WriteFile(path, []byte("a,b,c,d,count\n"+tc.rows), 0600); err != nil {
				t.Fatal(err)
			}
			counts, err := readQuartetCounts(path, tre)
			if (err != nil) != tc.fail {
				t.Fatalf("error=%v, expected failure=%v", err, tc.fail)
			}
			if tc.name == "valid" {
				var sum uint32
				for _, c := range counts {
					sum += c
				}
				if sum != 10 || len(counts) != 3 {
					t.Fatal(counts)
				}
			}
			if tc.name == "equivalent" {
				if len(counts) != 1 {
					t.Fatal(counts)
				}
				for _, c := range counts {
					if c != 7 {
						t.Fatal(counts)
					}
				}
			}
			if tc.name == "zeros" && len(counts) != 0 {
				t.Fatal(counts)
			}
		})
	}
}
