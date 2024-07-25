#include "si_calibration.h"

void create_gated_histograms(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    TH1D* histEnergySum[40][8];
    TH2D* histLeftRightGate[40][8][2];

    double gatingRange[40][8][fNumGates][2] = {0}; // det, strip, gate, range

    auto spectrum = new TSpectrum(fNumGates*2);
    TF1 *fitGaus = new TF1("fitGaus","gaus(0)",0,6000);

    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);

    for (auto det=0; det<40; ++det)
    {
        for (auto strip=0; strip<8; ++strip)
        {
            for (auto gate=0; gate<fNumGates; ++gate)
            {
                TString title = Form("detector %d, strip %d, gate %d", det, strip, gate);
                histLeftRightGate[det][strip][gate] = new TH2D(MakeHistName("LeftRightGated",det,strip,gate),title+";left;right",fNBinsA,fBinA1,fBinA2,fNBinsA,fBinA1,fBinA2);
                histLeftRightGate[det][strip][gate] -> SetStats(0);
            }
        }
    }

    auto fileHist = new TFile(fHistFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip) {
            histEnergySum[det][strip] = (TH1D*) fileHist -> Get(MakeHistName("EnergySum",det,strip));
        }

    for (auto det=0; det<40; ++det)
    {
        TCanvas *cvs = nullptr;

        bool drawCurrentDet = false;
        if (drawExample && (det==fExampleDet1||det==fExampleDet2))
            drawCurrentDet = true;

        if (drawCurrentDet) {
            cvs = MakeCanvas("Energy",det);
        }
        for (auto strip=0; strip<8; ++strip)
        {
            auto hist = histEnergySum[det][strip];
            if (drawCurrentDet) {
                cvs -> cd(strip+1);
                hist -> Draw();
            }
            auto numPeaks = spectrum -> Search(hist,5,"nodraw");
            double* xPeaks = spectrum -> GetPositionX();
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
            {
                if (iPeak+1>numPeaks) {
                    e_warning << "(" << det << ", " << strip << ") : " << " number of peaks is " << numPeaks << endl;
                    gatingRange[det][strip][iPeak][0] = 0;
                    gatingRange[det][strip][iPeak][1] = 0;
                    continue;
                }
                double xPeak = xPeaks[iPeak];
                fitGaus -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
                fitGaus -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
                hist -> Fit(fitGaus,"Q0NR");
                gatingRange[det][strip][iPeak][0] = fitGaus -> GetParameter(1) - 3 * fitGaus -> GetParameter(2);
                gatingRange[det][strip][iPeak][1] = fitGaus -> GetParameter(1) + 3 * fitGaus -> GetParameter(2);

                if (drawCurrentDet) {
                    cvs -> cd(strip+1);
                    fitGaus -> DrawClone("samel");
                }
            }
        }
    }

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0)
            cout << iEvent << " / " << numEvents << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto strip = channel -> GetStrip();
            auto energy1 = channel -> GetEnergy();
            auto energy2 = channel -> GetEnergy2();
            if (energy2>0) {
                auto sum = energy1 + energy2;
                auto pos = (energy1 - energy2) / sum;
                for (auto gate=0; gate<fNumGates; ++gate) {
                    if (sum>gatingRange[det][strip][gate][0] && sum<gatingRange[det][strip][gate][1])
                        histLeftRightGate[det][strip][gate] -> Fill(energy1, energy2);
                }
            }
        }
    }

    auto fileGate = new TFile(fGateFileName,"recreate");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip)
            for (auto gate=0; gate<fNumGates; ++gate)
                histLeftRightGate[det][strip][gate] -> Write();

    cout << fileGate -> GetName() << endl;
}
