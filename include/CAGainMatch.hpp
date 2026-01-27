#ifndef GAINMATCH_HPP
#define GAINMATCH_HPP

// C++ includes
#include <iostream>
#include <memory>
#include <vector>

// ROOT Includes
#include <TF1.h>
#include <TH2D.h>
#include <TSpectrum.h>

// Project Includes

namespace GainMatchConfig
{
    // Histogram Config
    inline const char* kAmplitudeHistogramName = "clover_cross/cc_amp"; // Path to the amplitude histogram within the root file
    inline constexpr unsigned int kNumChannels = 16; // Number of channels in the histogram

    // Background Subtraction Config
    inline constexpr unsigned int kBackgroundSmoothing = 50; // Smoothing parameter for background subtraction
    inline const char* kBackgroundOptions = "Compton"; // Sigma for peak finding

    // Peak Finding Config
    inline constexpr unsigned int kRebinFactor = 4; // Rebin factor for histograms
    inline constexpr unsigned int kMaxPeaks = 10; // Maximum number of peaks to identify
    inline constexpr std::pair<double, double> kPeakSearchRange = { 6000, 22000.0 }; // Bin range in which to search for peaks
    inline constexpr double kPeakSigma = 15.0; // Expected peak width (sigma) in bins
    inline constexpr double kPeakThreshold = 0.10; // Minimum height (as fraction of max) to consider a peak

    // Peak Matching Config
    inline constexpr std::pair<double, double> kReferenceEnergies = { 1460.820, 2614.511 }; // The background peaks we wish to use as our basis for gain matching
    inline constexpr double kPeakCentroidRatio = kReferenceEnergies.second / kReferenceEnergies.first; // Ratio describing 208Tl and 40K background lines
    inline constexpr double kPeakRatioTolerance = 0.0025; // Acceptable deviation in ratio when matching peaks

    // Fitting Config
    inline constexpr double kFitBounds = 3; // How many sigma to include in fit range

} // namespace GainMatchConfig

TH1D* BackgroundSubtraction1D(TH1D* hist)
{
    using namespace GainMatchConfig;
    auto spectrum = new TSpectrum(kMaxPeaks);
    auto bg =
        spectrum->Background(hist, kBackgroundSmoothing, kBackgroundOptions);
    hist->Add(bg, -1);

    delete bg;
    delete spectrum;

    return hist;
}

TH2D* BackgroundSubtraction2D(TH2D* hist)
{
    for (size_t ch = 0; ch < hist->GetNbinsY();
        ++ch) // Loop over each bin in Y (channel)
    {
        auto proj = hist->ProjectionX("_px", ch + 1, ch + 1);
        proj = BackgroundSubtraction1D(proj);

        for (int i = 0; i <= proj->GetNbinsX();
            ++i) // Copy back to original 2D histogram
            hist->SetBinContent(ch + 1, i, proj->GetBinContent(i));
        delete proj;
    }

    return hist;
}

std::vector<double> FindAllPeaks1D(TH1D* hist)
{
    using namespace GainMatchConfig;
    std::vector<double> all_peaks; // All peaks found in the search range

    auto spectrum = new TSpectrum(kMaxPeaks);

    hist->GetXaxis()->SetRangeUser(kPeakSearchRange.first,
        kPeakSearchRange.second); // Set search range
    auto n_found = spectrum->Search(hist, kPeakSigma, "", kPeakThreshold);
    auto peak_pos = spectrum->GetPositionX();

    all_peaks.reserve(n_found);
    for (size_t p = 0; p < n_found; ++p)
    {
        all_peaks.push_back(peak_pos[p]);
    }
    std::sort(all_peaks.begin(), all_peaks.end());

#if DEBUG >= 2
    printf("Found %u peak(s):\n", n_found); // Channel 0 is in bin 1
    for (size_t p = 0; p < all_peaks.size(); ++p)
    {
        printf("%.3f, ", all_peaks[p]); // For whatever reason we don't have to
        // multiply by kRebinFactor here
    }
    printf("\n");
#endif

    delete spectrum;

    return all_peaks;
}

std::vector<std::vector<double>> FindAllPeaks2D(TH2D* hist)
{
    std::vector<std::vector<double>> all_peaks_per_channel(
        GainMatchConfig::kNumChannels);

    for (size_t ch = 0; ch < hist->GetNbinsY(); ++ch) // Loop over each channel
    {
    #if DEBUG >= 2
        printf("Finding peaks in channel %zu (bin %zu)\n", ch,
            ch + 1); // Channel 0 is in bin 1
    #endif
        auto hist_proj = hist->ProjectionX("_px", ch + 1, ch + 1);
        auto channel_peaks = FindAllPeaks1D(hist_proj);
        all_peaks_per_channel[ch] = channel_peaks;
        delete hist_proj;
    }

    return all_peaks_per_channel;
}

std::pair<double, double> FindMatchingPeaks1D(std::vector<double>& peaks)
{
    std::pair<double, double> matching_peaks;

    double best_ratio = 0.0;
    auto best_pair = std::make_pair<double, double>(-1.0, -1.0);

    // Look through all combinations of peaks for this channel
    using namespace GainMatchConfig;
    for (size_t i = 0; i < peaks.size(); ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            auto ratio = peaks[i] / peaks[j];
            if (std::abs(ratio - kPeakCentroidRatio) / kPeakCentroidRatio <
                kPeakRatioTolerance) // Allow tolerance about the expected ratio
            {
                if (std::abs(ratio - kPeakCentroidRatio) <
                    std::abs(best_ratio - kPeakCentroidRatio))
                {
                    best_ratio = ratio;
                    best_pair = std::make_pair(
                        peaks[j], peaks[i]); // Store as (low, high)
                #if DEBUG >= 2
                    printf("Found matching peaks: (%.0f,%.0f) with ratio match "
                        "(%f)\n",
                        peaks[i], peaks[j], ratio / kPeakCentroidRatio);
                #endif
                }
            }
        }
    }
    if (best_pair.first > 0 && best_pair.second > 0)
    {
        matching_peaks = best_pair;
    }
    else
    {
        printf("No matching peaks found within tolerance (%f)\n",
            kPeakRatioTolerance);
    }

    return matching_peaks;
}

std::vector<std::pair<double, double>>
FindMatchingPeaks2D(std::vector<std::vector<double>>& peaks)
{
    std::vector<std::pair<double, double>> all_matching_peaks_per_channel(
        GainMatchConfig::kNumChannels);
    for (size_t ch = 0; ch < peaks.size(); ++ch)
    {
    #if DEBUG >= 2
        printf("Finding matching peaks in channel %zu\n", ch);
    #endif
        auto channel_matching_peaks = FindMatchingPeaks1D(peaks[ch]);
        all_matching_peaks_per_channel[ch] = channel_matching_peaks;
    }
    return all_matching_peaks_per_channel;
}

std::pair<double, double>
GetPeakCentroids1D(TH1D* hist, std::pair<double, double>& matched_peaks)
{
    using namespace GainMatchConfig;

    std::pair<double, double> centroid_pair(-1.0, -1.0);
    bool got_first = false;

    for (auto& peak_pos : { matched_peaks.first, matched_peaks.second })
    {
        double peak_height = hist->GetBinContent(hist->FindBin(peak_pos));
        double fit_min = peak_pos - kFitBounds * kPeakSigma;
        double fit_max = peak_pos + kFitBounds * kPeakSigma;

        auto gaus = new TF1("gaus", "gaus", fit_min, fit_max);

        // Initial Values for Fit
        gaus->SetParameters(peak_height, peak_pos, kPeakSigma);
        gaus->SetParLimits(1, peak_pos - kPeakSigma, peak_pos + kPeakSigma);
        gaus->SetParLimits(2, 0, kPeakSigma * 2.0);

        hist->Fit(gaus, "QLMRES0", "", fit_min, fit_max);
        auto mean = gaus->GetParameter(1);

        if (!got_first)
        {
            centroid_pair.first = mean;
            got_first = true;
        }
        else
        {
            centroid_pair.second = mean;
        }

    #if (DEBUG >= 2)
        printf("Fitted peak at %.3f with centroid %.3f (initial pos %.3f)\n",
            peak_pos, mean, peak_pos);
    #endif
        delete gaus;
    }

    return centroid_pair;
}

std::vector<std::pair<double, double>>
GetPeakCentroids2D(TH2D* hist,
    std::vector<std::pair<double, double>>& matched_peaks)
{
    std::vector<std::pair<double, double>> all_centroids_per_channel(
        GainMatchConfig::kNumChannels);

    for (size_t ch = 0; ch < matched_peaks.size(); ++ch)
    {
    #if DEBUG >= 2
        printf("Getting centroids for channel %zu\n", ch);
    #endif
        auto hist_proj = hist->ProjectionX("_px", ch + 1, ch + 1);
        auto channel_centroids =
            GetPeakCentroids1D(hist_proj, matched_peaks[ch]);
        all_centroids_per_channel[ch] = channel_centroids;
        delete hist_proj;
    }
    return all_centroids_per_channel;
}

std::vector<std::pair<double, double>> CalculateGainMatchParameters(
    std::vector<std::pair<double, double>>& ref_centroids,
    std::vector<std::pair<double, double>>& inp_centroids)
{
    std::vector<std::pair<double, double>> all_gainmatch_params; // has (gain, offset) pair for every channel

    if (ref_centroids.size() != inp_centroids.size())
    {
        std::cerr
            << "Error: Mismatched number of fitted peaks between reference "
            "and input!"
            << std::endl;
        return all_gainmatch_params;
    }

    size_t n_peaks = std::min(ref_centroids.size(), inp_centroids.size());
    for (size_t i = 0; i < n_peaks; ++i)
    {
        double ref_low = ref_centroids[i].first;
        double ref_high = ref_centroids[i].second;
        double input_low = inp_centroids[i].first;
        double input_high = inp_centroids[i].second;
        double gain = (ref_high - ref_low) / (input_high - input_low);
        double offset = ref_low - gain * input_low;
        all_gainmatch_params.push_back(std::make_pair(offset, gain));

        printf("Channel %zu: Offset = %.6f, Gain = %.10f\n", i, offset, gain);
    }
    return all_gainmatch_params;
}

void WriteParamatersToFile(
    const std::string& output_filename,
    const std::vector<std::pair<double, double>>& params)
{
    FILE* out_file = fopen(output_filename.c_str(), "w");
    if (!out_file)
    {
        std::cerr << "Error opening output file for writing" << std::endl;
        return;
    }

    fprintf(out_file, "# Channel\tOffset\tGain\n");
    fprintf(out_file, "# Clover Cross\n");
    for (size_t ch = 0; ch < params.size(); ++ch)
    {
        fprintf(out_file, "%zu\t%.10f\t%.10f\n", ch, params[ch].first,
            params[ch].second);
    }
    fclose(out_file);
    printf("Gain match parameters written to %s!\n", output_filename.c_str());
}

#endif // GAINMATCH_HPP