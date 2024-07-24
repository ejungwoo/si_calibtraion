#include "calibration.h"

void create_histograms()
{
    MakeRun(253);

    TH1D* histEnergySum[40][8];
    TH2D* histLeftRight[40][8];
    TH2D* histEnergyPosition[40][8];

    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);

    for (auto det=0; det<40; ++det)
    {
        for (auto strip=0; strip<8; ++strip)
        {
            TString title = Form("detector %d, strip %d", det, strip);
            histEnergySum[det][strip] = new TH1D(MakeHistName("EnergySum",det,strip),title+";energy;",fNBinsE,fBinE1,fBinE2);
            histEnergySum[det][strip] -> SetStats(0);
            histLeftRight[det][strip] = new TH2D(MakeHistName("LeftRight",det,strip),title+";left;right",fNBinsE,fBinE1,fBinE2,fNBinsE,fBinE1,fBinE2);
            histLeftRight[det][strip] -> SetStats(0);
            histEnergyPosition[det][strip] = new TH2D(MakeHistName("EnergyPosition",det,strip),title+";position;energy sum",200,-1,1,fNBinsE,fBinE1,1.5*fBinE2);
            histEnergyPosition[det][strip] -> SetStats(0);
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
                histEnergySum[det][strip] -> Fill(sum);
                histLeftRight[det][strip] -> Fill(energy1, energy2);
                histEnergyPosition[det][strip] -> Fill(pos,sum);
            }
        }
    }

    auto fileHist = new TFile(fHistFileName,"recreate");
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) histEnergySum[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) histLeftRight[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) histEnergyPosition[det][strip] -> Write();

    cout << fileHist -> GetName() << endl;
}
