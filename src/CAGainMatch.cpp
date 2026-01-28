// C++ includes
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

// ROOT Includes
#include <TF1.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TSpectrum.h>

// Project Includes
#include "CAGainMatch.hpp"

int main(int argc, char* argv[])
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

    /* #region Open Files */

    TFile* ref_file = nullptr;
    ReferenceFileType reference_file_type;
    if (reference_filename.rfind(".root") != std::string::npos)
    {
        printf("Opening reference ROOT file: %s\n", reference_filename.c_str());
        ref_file = TFile::Open(reference_filename.c_str(), "READ");
        if (!ref_file)
        {
            std::cerr << "Error opening reference file" << std::endl;
            return 1;
        }
        reference_file_type = ReferenceFileType::ROOT;
    }
    else if (reference_filename.rfind(".capk") != std::string::npos) // Clover Array Peaks File
    {
        printf("Reference file is a Clover Array Peak file, will load peaks from it when needed\n");
        reference_file_type = ReferenceFileType::CAPK;
    }
    else
    {
        std::cerr << "Error: Unsupported reference file format. Please provide a \".root\" or \".capk\" file." << std::endl;
        return 1;
    }

    TFile* inp_file = nullptr;
    if (input_filename.rfind(".root") != std::string::npos)
    {
        inp_file = TFile::Open(input_filename.c_str(), "READ");
        if (!inp_file)
        {
            std::cerr << "Error opening input file" << std::endl;
            return 1;
        }
    }
    else if (input_filename.rfind(".cags") != std::string::npos)
    {
        std::cerr << "Error: .cags files are not supported as input files. Please provide a .root file." << std::endl;
        return 1;
    }
    else
    {
        std::cerr << "Error: Unsupported input file format. Please provide a .root file." << std::endl;
        return 1;
    }

    /* #endregion Open Files */

    /* #region Get Histograms */
    TH2D* ref_hist = nullptr;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        ref_hist = dynamic_cast<TH2D*>(ref_file->Get(kAmplitudeHistogramName));
        if (!ref_hist)
        {
            std::cerr << "Error retrieving reference histogram" << std::endl;
            return 1;
        }
        ref_hist->RebinX(kRebinFactor); // Rebin to make peakfinding easier
    }

    auto inp_hist = dynamic_cast<TH2D*>(inp_file->Get(kAmplitudeHistogramName));
    if (!inp_hist)

    {
        std::cerr << "Error retrieving input histogram" << std::endl;
        return 1;
    }
    inp_hist->RebinX(kRebinFactor); // Rebin to make peakfinding easier

    /* #endregion Get Histograms */

    // Perform background subtraction to give peak finding a better chance
    printf("Performing background subtraction...\n");
    if (reference_file_type == ReferenceFileType::ROOT)
        BackgroundSubtraction2D(ref_hist);
    BackgroundSubtraction2D(inp_hist);

    // Find Peaks
    std::vector<std::vector<double>> ref_all_peaks;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        printf("Finding peaks in reference histogram: %s\n", ref_hist->GetName());
        ref_all_peaks = FindAllPeaks2D(ref_hist);
    }
    printf("Finding peaks in input histogram: %s\n", inp_hist->GetName());
    auto inp_all_peaks = FindAllPeaks2D(inp_hist);

    // Find matching peaks based on expected centroid ratio
    std::vector<std::pair<double, double>> ref_matched_peaks;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        printf("Finding matching peaks in reference histogram...\n");
        ref_matched_peaks = FindMatchingPeaks2D(ref_all_peaks);
    }
    printf("Finding matching peaks in input histogram...\n");
    auto inp_matched_peaks = FindMatchingPeaks2D(inp_all_peaks);

    // Fit the matched peaks to get more precise centroids
    std::vector<std::pair<double, double>> ref_centroids;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        printf("Fitting matched peaks in reference histogram...\n");
        ref_centroids = GetPeakCentroids2D(ref_hist, ref_matched_peaks);
    }
    else if (reference_file_type == ReferenceFileType::CAPK)
    {
        printf("Loading peak centroids from CAPK file: %s\n", reference_filename.c_str());
        ref_centroids = LoadPeaksFromCAPKFile(reference_filename);
    }


    printf("Fitting matched peaks in input histogram...\n");
    auto inp_centroids = GetPeakCentroids2D(inp_hist, inp_matched_peaks);

    printf("Obtained %zu reference centroid and %zu input centroid pairs.\n",
        ref_centroids.size(), inp_centroids.size());
    printf("Centroids (Reference, Input): \n");
    for (size_t ch = 0; ch < inp_centroids.size(); ++ch)
    {
        printf("Channel %zu: (%.3f, %.3f) , (%.3f, %.3f)\n", ch,
            ref_centroids[ch].first, ref_centroids[ch].second,
            inp_centroids[ch].first, inp_centroids[ch].second);

    }

    // Calculate Gain Match Parameters
    printf("Calculating gain match parameters...\n");
    auto params = CalculateGainMatchParameters(ref_centroids, inp_centroids);

    // Output Gain Match Parameters to Output File
    std::string final_output_filename = output_filename;
    if (output_filename.size() < 5 || output_filename.substr(output_filename.size() - 5) != ".cags")
    {
        size_t dot_pos = output_filename.rfind('.');
        if (dot_pos != std::string::npos)
        {
            final_output_filename = output_filename.substr(0, dot_pos) + ".cags";
        }
        else
        {
            final_output_filename = output_filename + ".cags";
        }
    }

    FILE* out_file = fopen(final_output_filename.c_str(), "w");
    if (!out_file)
    {
        std::cerr << "Error opening output file for writing" << std::endl;
        return 1;
    }

    fprintf(out_file, "# Channel\tOffset\tGain\n");
    fprintf(out_file, "# Clover Cross\n");
    for (size_t ch = 0; ch < params.size(); ++ch)
    {
        fprintf(out_file, "%zu\t%.10f\t%.10f\n", ch, params[ch].first,
            params[ch].second);
    }
    fclose(out_file);
    printf("Gain match parameters written to %s!\n", final_output_filename.c_str());

    // Close input files
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        ref_file->Close();
        delete ref_file;
    }
    inp_file->Close();
    delete inp_file;

    printf("Gain matching completed successfully!\n");

    return 0;
}
