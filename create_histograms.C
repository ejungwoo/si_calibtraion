#include "si_calibration.h"

void create_histograms(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);

    for (auto det=0; det<40; ++det)
    {
        for (auto strip=0; strip<8; ++strip)
        {
            TString title = Form("detector %d, strip %d", det, strip);
            fHistEnergySum[det][strip] = new TH1D(MakeHistName("EnergySum",det,strip),title+";energy;",fNBinsA,fBinA1,fBinA2);
            fHistEnergySum[det][strip] -> SetStats(0);
            fHistLeftRight[det][strip] = new TH2D(MakeHistName("LeftRight",det,strip),title+";left;right",fNBinsA,fBinA1,fBinA2,fNBinsA,fBinA1,fBinA2);
            fHistLeftRight[det][strip] -> SetStats(0);
            fHistEnergyPosition[det][strip] = new TH2D(MakeHistName("EnergyPosition",det,strip),title+";position;energy sum",200,-1,1,fNBinsA,fBinA1,1.5*fBinA2);
            fHistEnergyPosition[det][strip] -> SetStats(0);
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
                fHistEnergySum[det][strip] -> Fill(sum);
                fHistLeftRight[det][strip] -> Fill(energy1, energy2);
                fHistEnergyPosition[det][strip] -> Fill(pos,sum);
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
            auto cvs1 = MakeCanvas("EnergySum",det);
            auto cvs2 = MakeCanvas("LeftRight",det);
            auto cvs3 = MakeCanvas("EnergyPosition",det);
            for (auto strip=0; strip<8; ++strip) {
                cvs1 -> cd(strip+1); fHistEnergySum[det][strip] -> Draw();
                cvs2 -> cd(strip+1); fHistLeftRight[det][strip] -> Draw("colz");
                cvs3 -> cd(strip+1); fHistEnergyPosition[det][strip] -> Draw("colz");
            }
        }
    }

    auto fileHist = new TFile(fHistFileName,"recreate");
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistEnergySum[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistLeftRight[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistEnergyPosition[det][strip] -> Write();

    cout << fileHist -> GetName() << endl;
}
