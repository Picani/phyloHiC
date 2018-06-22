package main

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"

	"gonum.org/v1/gonum/stat"
	"math"
	"sort"
)

func ReadValuesNoZeros(filename string) ([]float64, error) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		return nil, err
	}

	var reader *csv.Reader
	if filepath.Ext(filename) == ".gz" {
		comp, err := gzip.NewReader(f)
		defer comp.Close()
		if err != nil {
			return nil, err
		}
		reader = csv.NewReader(comp)
	} else {
		reader = csv.NewReader(f)
	}
	reader.Comma = '\t'

	vals := make([]float64, 0)

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}

		x, err := strconv.ParseFloat(record[2], 64)
		if err != nil {
			return nil, err
		}
		if x == 0.0 {
			continue
		}

		vals = append(vals, x)
	}

	return vals, nil
}

func Median(x []float64) float64 {
	if !sort.Float64sAreSorted(x) {
		return math.NaN()
	}
	if len(x) == 0 {
		return math.NaN()
	}

	if len(x)%2 != 0 {
		return x[len(x)/2]
	} else {
		a := x[(len(x)/2)-1]
		b := x[len(x)/2]
		fmt.Printf("a: %v\tb: %v\n", a, b)
		return (a + b) / 2.0
	}
}

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "No! Please God No!")
		os.Exit(1)
	}

	means := make(map[string][]float64)
	means["Size"] = make([]float64, len(os.Args)-1)
	means["Mean"] = make([]float64, len(os.Args)-1)
	means["StdErr"] = make([]float64, len(os.Args)-1)
	means["StdDev"] = make([]float64, len(os.Args)-1)
	means["10%"] = make([]float64, len(os.Args)-1)
	means["25%"] = make([]float64, len(os.Args)-1)
	means["Median"] = make([]float64, len(os.Args)-1)
	means["75%"] = make([]float64, len(os.Args)-1)
	means["90%"] = make([]float64, len(os.Args)-1)

	fmt.Println("Name\tSize\tMean\tStdErr\tStdDev\t10%\t25%\tMedian\t75%\t90%")
	for i, file := range os.Args[1:] {
		vals, err := ReadValuesNoZeros(file)
		if err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		sort.Float64s(vals)
		mean := stat.Mean(vals, nil)
		stddev := stat.StdDev(vals, nil)
		stderr := stat.StdErr(stddev, float64(len(vals)))
		q10 := stat.Quantile(0.1, 1, vals, nil)
		q25 := stat.Quantile(0.25, 1, vals, nil)
		median := stat.Quantile(0.5, 1, vals, nil)
		q75 := stat.Quantile(0.75, 1, vals, nil)
		q90 := stat.Quantile(0.9, 1, vals, nil)
		fmt.Printf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", file, len(vals), mean, stderr, stddev, q10, q25, median, q75, q90)

		means["Size"][i] = float64(len(vals))
		means["Mean"][i] = mean
		means["StdErr"][i] = stderr
		means["StdDev"][i] = stddev
		means["10%"][i] = q10
		means["25%"][i] = q25
		means["Median"][i] = median
		means["75%"][i] = q75
		means["90%"][i] = q90
	}

	fmt.Printf("%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		"Means",
		stat.Mean(means["Size"], nil),
		stat.Mean(means["Mean"], nil),
		stat.Mean(means["StdErr"], nil),
		stat.Mean(means["StdDev"], nil),
		stat.Mean(means["10%"], nil),
		stat.Mean(means["25%"], nil),
		stat.Mean(means["Median"], nil),
		stat.Mean(means["75%"], nil),
		stat.Mean(means["90%"], nil))

	os.Exit(0)
}
