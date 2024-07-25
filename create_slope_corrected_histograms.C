#include "si_calibration.h"

void create_slope_corrected_histograms(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    ifstream fcalibration(fSlopeCalibrationParFileName);
    int det, strip;
    double g1, g2;
    double g1Array[40][8];
    double g2Array[40][8];
    while (fcalibration >> det >> strip >> g1 >> g2)
    {
        g1Array[det][strip] = abs(g1);
        g2Array[det][strip] = abs(g2);
    }

    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);

    for (auto det=0; det<40; ++det)
    {
        for (auto strip=0; strip<8; ++strip)
        {
            TString title = Form("detector %d, strip %d", det, strip);
            fHistLeftRightSlopeCalibrated[det][strip] = new TH2D(MakeHistName("LeftRightSlopeCalibrated",det,strip),title+";left;right",fNBinsE,fBinE1,fBinE2,fNBinsE,fBinE1,fBinE2);
            fHistLeftRightSlopeCalibrated[det][strip] -> SetStats(0);
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
            auto energy1 = channel -> GetEnergy() * g1Array[det][strip];
            auto energy2 = channel -> GetEnergy2() * g2Array[det][strip];
            if (energy2>0) {
                auto sum = energy1 + energy2;
                auto pos = (energy1 - energy2) / sum;
                fHistLeftRightSlopeCalibrated[det][strip] -> Fill(energy1, energy2);
            }
        }
    }

    for (auto det=0; det<40; ++det)
    {
        TCanvas *cvs = nullptr;

        bool drawCurrentDet = false;
        if (drawExample && (det==fExampleDet1||det==fExampleDet2))
            drawCurrentDet = true;

        if (drawCurrentDet) {
            auto cvs = MakeCanvas("LeftRightSlopeCalibrated",det);
            for (auto strip=0; strip<8; ++strip) {
                cvs -> cd(strip+1);
                fHistLeftRightSlopeCalibrated[det][strip] -> Draw("colz");
            }
        }
    }

    auto fileHist = new TFile(fSlopeCalibratedHistFileName,"recreate");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip)
            fHistLeftRightSlopeCalibrated[det][strip] -> Write();

    cout << fileHist -> GetName() << endl;
}
