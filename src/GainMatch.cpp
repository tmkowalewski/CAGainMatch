// C++ includes
#include <iostream>
#include <memory>

// ROOT Includes
#include <TFile.h>
#include <TH2D.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TF1.h>

// Project Includes
#include "GainMatch.hpp"

// Configuration
#define DEBUG 0

int main(int argc, char *argv[])
{
    // Introduction
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <reference file> <input file> <output file>" << std::endl;
        return 1;
    }

    const std::string reference_filename = argv[1];
    const std::string input_filename = argv[2];
    const std::string output_filename = argv[3];
    printf("============= Welcome to GainMatch! ===============\n");
    printf("------------- Current Configuration ---------------\n");
    printf("Using reference file: %s\n", reference_filename.c_str());
    printf("Using input file: %s\n", input_filename.c_str());
    printf("Output file: %s\n", output_filename.c_str());
    printf("---------------------------------------------------\n");

    std::cout << "Gain matching started!" << std::endl;

    using namespace GainMatchConfig;

    // Open Files

    auto ref_file = TFile::Open(reference_filename.c_str(), "READ");
    if (!ref_file)
    {
        std::cerr << "Error opening reference file" << std::endl;
        return 1;
    }

    auto inp_file = TFile::Open(input_filename.c_str(), "READ");
    if (!inp_file)
    {
        std::cerr << "Error opening input file" << std::endl;
        return 1;
    }

    // Get Histograms

    auto ref_hist = dynamic_cast<TH2D *>(ref_file->Get(kAmplitudeHistogramName));
    if (!ref_hist)
    {
        std::cerr << "Error retrieving reference histogram" << std::endl;
        return 1;
    }
    ref_hist->RebinX(kRebinFactor); // Rebin to make peakfinding easier

    auto inp_hist = dynamic_cast<TH2D *>(inp_file->Get(kAmplitudeHistogramName));
    if (!inp_hist)

    {
        std::cerr << "Error retrieving input histogram" << std::endl;
        return 1;
    }
    inp_hist->RebinX(kRebinFactor); // Rebin to make peakfinding easier

    // Perform background subtraction to give peak finding a better chance
    printf("Performing background subtraction...\n");
    BackgroundSubtraction2D(ref_hist);
    BackgroundSubtraction2D(inp_hist);

    // Find Peaks
    printf("Finding peaks in reference histogram: %s\n", ref_hist->GetName());
    auto ref_all_peaks = FindAllPeaks2D(ref_hist);
    printf("Finding peaks in input histogram: %s\n", inp_hist->GetName());
    auto inp_all_peaks = FindAllPeaks2D(inp_hist);

    // Find matching peaks based on expected centroid ratio
    printf("Finding matching peaks in reference histogram...\n");
    auto ref_matched_peaks = FindMatchingPeaks2D(ref_all_peaks);
    printf("Finding matching peaks in input histogram...\n");
    auto inp_matched_peaks = FindMatchingPeaks2D(inp_all_peaks);

    // Fit the matched peaks to get more precise centroids
    printf("Fitting matched peaks in reference histogram...\n");
    auto ref_centroids = GetPeakCentroids2D(ref_hist, ref_matched_peaks);
    printf("Fitting matched peaks in input histogram...\n");
    auto inp_centroids = GetPeakCentroids2D(inp_hist, inp_matched_peaks);

    // Calculate Gain Match Parameters
    printf("Calculating gain match parameters...\n");
    auto params = CalculateGainMatchParameters(ref_centroids, inp_centroids);

    // Close input files
    ref_file->Close();
    delete ref_file;
    inp_file->Close();
    delete inp_file;

    return 0;
}
