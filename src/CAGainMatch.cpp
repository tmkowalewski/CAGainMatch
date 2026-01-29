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

    printf("[INFO] Gain matching started!\n");

    using namespace GainMatchConfig;

    /* #region Open Files */

    TFile* ref_file = nullptr;
    ReferenceFileType reference_file_type;
    if (reference_filename.rfind(".root") != std::string::npos)
    {
        printf("[INFO] Opening reference ROOT file: %s\n", reference_filename.c_str());
        ref_file = TFile::Open(reference_filename.c_str(), "READ");
        if (!ref_file)
        {
            throw std::runtime_error(Form("[ERROR] Could not open reference ROOT file: %s", reference_filename.c_str()));
        }
        reference_file_type = ReferenceFileType::ROOT;
    }
    else if (reference_filename.rfind(".capk") != std::string::npos) // Clover Array Peaks File
    {
        printf("[INFO] Reference file is a Clover Array Peak file, will load peaks from it when needed\n");
        reference_file_type = ReferenceFileType::CAPK;
    }
    else
    {
        throw std::runtime_error(Form("[ERROR] Unsupported reference file format: %s", reference_filename.c_str()));
    }

    TFile* inp_file = nullptr;
    if (input_filename.rfind(".root") != std::string::npos)
    {
        inp_file = TFile::Open(input_filename.c_str(), "READ");
        if (!inp_file)
        {
            throw std::runtime_error(Form("[ERROR] Could not open input ROOT file: %s", input_filename.c_str()));
        }
    }
    else
    {
        throw std::runtime_error("[ERROR] Unsupported input file format. Please provide a .root file.");
    }

    /* #endregion Open Files */

    /* #region Get Histograms */
    std::vector<TH2D*> ref_hists;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        for (auto& hist_name : GainMatchConfig::kAmplitudeHistogramNames)
        {
            ref_hists.push_back(dynamic_cast<TH2D*>(ref_file->Get(hist_name.c_str())));
            if (!ref_hists.back())
            {
                throw std::runtime_error(Form("[ERROR] Could not retrieve reference histogram: %s\n", hist_name.c_str()));
            }
            ref_hists.back()->RebinX(kRebinFactor); // Rebin to make peakfinding easier
        }
    }

    std::vector<TH2D*> inp_hists;
    for (auto& hist_name : GainMatchConfig::kAmplitudeHistogramNames)
    {
        inp_hists.push_back(dynamic_cast<TH2D*>(inp_file->Get(hist_name.c_str())));
        if (!inp_hists.back())
        {
            throw std::runtime_error(Form("Error retrieving input histogram: %s\n", hist_name.c_str()));
        }
        inp_hists.back()->RebinX(kRebinFactor); // Rebin to make peakfinding easier
    }

    /* #endregion Get Histograms */

    /* #region Background Subtraction */
    printf("[INFO] Performing background subtraction\n");
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        for (auto& hist : ref_hists)
        {
            BackgroundSubtraction2D(hist);
        }
    }

    for (auto& hist : inp_hists)
        BackgroundSubtraction2D(hist);

    /* #endregion Background Subtraction */

    /* #region Peak Finding */
    std::vector<std::vector<std::vector<double>>> ref_all_peaks, inp_all_peaks;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        for (auto& hist : ref_hists)
        {
            printf("[INFO] Finding peaks in reference histogram: %s\n", hist->GetName());
            ref_all_peaks.push_back(FindAllPeaks2D(hist));
        }
    }

    for (auto& hist : inp_hists)
    {
        printf("[INFO] Finding peaks in input histogram: %s\n", hist->GetName());
        inp_all_peaks.push_back(FindAllPeaks2D(hist));
    }
    /* #endregion Peak Finding */

    /* #region Find Matching Peaks */
    std::vector<std::vector<std::pair<double, double>>> ref_matched_peaks, inp_matched_peaks;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        for (size_t i = 0; i < ref_all_peaks.size(); ++i)
        {
            printf("[INFO] Finding matching peaks in reference histogram %s\n", ref_hists[i]->GetName());
            ref_matched_peaks.push_back(FindMatchingPeaks2D(ref_all_peaks[i]));
        }
    }

    printf("[INFO] Finding matching peaks in input histogram...\n");
    for (size_t i = 0; i < inp_all_peaks.size(); ++i)
    {
        printf("[INFO] Finding matching peaks in input histogram %s\n", inp_hists[i]->GetName());
        inp_matched_peaks.push_back(FindMatchingPeaks2D(inp_all_peaks[i]));
    }
    /* #endregion Find Matching Peaks */

    /* #region Fit precise Peak Centroids */
    std::vector<std::vector<std::pair<double, double>>> ref_centroids, inp_centroids;
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        for (size_t i = 0; i < ref_matched_peaks.size(); ++i)
        {
            printf("[INFO] Fitting matched peaks in reference histogram %s...\n", ref_hists[i]->GetName());
            ref_centroids.push_back(GetPeakCentroids2D(ref_hists[i], ref_matched_peaks[i]));
        }
    }
    else if (reference_file_type == ReferenceFileType::CAPK)
    {
        printf("[INFO] Loading peak centroids from CAPK file: %s\n", reference_filename.c_str());
        ref_centroids = LoadPeaksFromCAPKFile(reference_filename);
    }

    printf("[INFO] Fitting matched peaks in input histogram...\n");
    for (size_t i = 0; i < inp_matched_peaks.size(); ++i)
    {
        printf("[INFO] Fitting matched peaks in input histogram %s...\n", inp_hists[i]->GetName());
        inp_centroids.push_back(GetPeakCentroids2D(inp_hists[i], inp_matched_peaks[i]));
    }
    for (size_t ch = 0; ch < kAmplitudeHistogramNames.size(); ++ch)
    {
        #if DEBUG >= 2
        printf("[DEBUG] Obtained %zu reference centroid and %zu input centroid pairs for histogram %s\n:",
            ref_centroids.size(), inp_centroids.size(), kAmplitudeHistogramNames[ch].c_str());
        printf("[DEBUG] Channel\tReference Pair\tInput Pair\n");
        for (size_t c = 0; c < ref_centroids.size(); ++c)
        {
            printf("[DEBUG] %zu\t(%.3f, %.3f)\t(%.3f, %.3f)\n", c,
                ref_centroids[c].first, ref_centroids[c].second,
                inp_centroids[c].first, inp_centroids[c].second);
        }
        #endif
    }
    /* #endregion Fit precise Peak Centroids */

    /* #region Calculate Gain Match Parameters */
    printf("[INFO] Calculating gain match parameters...\n");
    std::vector<std::vector<std::pair<double, double>>> all_params;
    for (size_t i = 0; i < kAmplitudeHistogramNames.size(); ++i)
    {
        printf("[INFO] Calculating gain match parameters for histogram %s\n", kAmplitudeHistogramNames[i].c_str());
        all_params.push_back(CalculateGainMatchParameters(ref_centroids[i], inp_centroids[i]));
    }
    /* #endregion Calculate Gain Match Parameters */

    /* #region Output Gain Match Parameters to Output File */
    std::string final_output_filename = output_filename;
    // Ensure output filename ends with .cags
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
        throw std::runtime_error(Form("[ERROR] Could not open output file for writing: %s", final_output_filename.c_str()));
    }

    fprintf(out_file, "# Channel\tOffset\tGain\n");
    for (size_t md = 0; md < all_params.size(); ++md)
    {
        std::string hist_name = kAmplitudeHistogramNames[md];
        size_t slash_pos = hist_name.find('/');
        if (hist_name.find('/') != std::string::npos)
        {
            hist_name = hist_name.substr(0, slash_pos);
        }
        fprintf(out_file, "# %s\n", hist_name.c_str());
        const auto& params = all_params[md];
        for (size_t ch = 0; ch < params.size(); ++ch)
        {
            fprintf(out_file, "%2d\t% 14.10f\t% 14.10f\n", ch, params[ch].first,
                params[ch].second);
        }
    }

    fclose(out_file);
    printf("[INFO] Gain match parameters written to %s!\n", final_output_filename.c_str());
    /* #endregion Output Gain Match Parameters to Output File */

    // Clean-up
    if (reference_file_type == ReferenceFileType::ROOT)
    {
        ref_file->Close();
        delete ref_file;
    }
    inp_file->Close();
    delete inp_file;

    printf("[INFO] Gain matching completed successfully!\n");

    return 0;
}
