package sgp4

import (
	"fmt"
	"os"
	"strings"
	"testing"
)

func TestOMMToTLEStringSamples(t *testing.T) {
	jsonData, err := os.ReadFile("test_samples/omms.json")
	if err != nil {
		t.Fatal(err)
	}
	messages, err := ParseOMMs(jsonData)
	if err != nil {
		t.Fatal(err)
	}
	tleData, err := os.ReadFile("test_samples/tle.txt")
	if err != nil {
		t.Fatal(err)
	}
	fixture := strings.TrimSuffix(string(tleData), "\n")
	lines := strings.Split(fixture, "\n")
	if len(lines) != 3*len(messages) {
		t.Fatalf("sample count mismatch: %d OMM records and %d TLE lines", len(messages), len(lines))
	}

	for i := range messages {
		message := &messages[i]
		t.Run(fmt.Sprintf("%02d_NORAD-%d", i+1, message.NoradCatID), func(t *testing.T) {
			tle, err := message.ToTLE()
			if err != nil {
				t.Fatalf("issOMM.ToTLE() failed: %v", err)
			}

			actual, err := tle.TLEString()
			if err != nil {
				t.Fatal(err)
			}
			actualLines := strings.Split(actual, "\n")
			if len(actualLines) != 3 {
				t.Fatalf("generated %d lines, want 3", len(actualLines))
			}
			if got, want := strings.TrimSpace(actualLines[0]), strings.TrimSpace(lines[3*i]); got != want {
				t.Errorf("name = %q, want %q", got, want)
			}
			for line := 1; line <= 2; line++ {
				if got, want := actualLines[line], lines[3*i+line]; got != want {
					t.Errorf("TLE line %d mismatch:\n got: %q\nwant: %q", line, got, want)
				}
			}
		})
	}
}
