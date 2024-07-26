#include "si_calibration.h"

void make_c2_calibration_histograms(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    int ne = fNBinsE;
    double e1 = fBinE1;
    double e2 = fBinE2;

    double gArray[40][8][2];
    double bArray[40][8][3];
    TString calibrationName = "C2";

    ifstream fcalibration_g(fC1ParFileName);
    int det, strip;
    double g1, g2;
    while (fcalibration_g >> det >> strip >> g1 >> g2)
    {
        gArray[det][strip][0] = abs(g1);
        gArray[det][strip][1] = abs(g2);
    }
    fcalibration_g.close();

    ifstream fcalibration_b(fC2ParFileName);
    double b0, b1, b2;
    while (fcalibration_b >> det >> strip >> b0 >> b1 >> b2)
    {
        bArray[det][strip][0] = b0;
        bArray[det][strip][1] = b1;
        bArray[det][strip][2] = b2;
    }
    fcalibration_b.close();

    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);

    fHistEnergyPositionAll = new TH2D(MakeHistName(calibrationName+"EnergyPositionAll"),";position;energy sum",400,-1,1,2*ne,e1,1.5*e2); 
    for (auto det=0; det<40; ++det)
    {
        for (auto strip=0; strip<8; ++strip)
        {
            TString title = Form("detector %d, strip %d", det, strip);
            fHistEnergySum[det][strip] = new TH1D(MakeHistName(calibrationName+"EnergySum",det,strip),title+";energy;",ne,e1,e2);
            fHistEnergySum[det][strip] -> SetStats(0);
            fHistLeftRight[det][strip] = new TH2D(MakeHistName(calibrationName+"LeftRight",det,strip),title+";left;right",ne,e1,e2,ne,e1,e2);
            fHistLeftRight[det][strip] -> SetStats(0);
            fHistEnergyPosition[det][strip] = new TH2D(MakeHistName(calibrationName+"EnergyPosition",det,strip),title+";position;energy sum",200,-1,1,ne,e1,1.5*e2);
            fHistEnergyPosition[det][strip] -> SetStats(0);
        }
    }

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0)
            cout << "Filling histogram " << iEvent << " / " << numEvents << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto strip = channel -> GetStrip();
            double g1 = gArray[det][strip][0];
            double g2 = gArray[det][strip][1];
            double b0 = bArray[det][strip][0];
            double b1 = bArray[det][strip][1];
            double b2 = bArray[det][strip][2];
            auto energy1 = channel -> GetEnergy() * g1;
            auto energy2 = channel -> GetEnergy2() * g2;
            if (energy2>0) {
                auto sum = (energy1 + energy2);
                auto pos = (energy1 - energy2) / sum;
                sum = sum / (b0 + b1*pos + b2*pos*pos) * fAlphaEnergy241Am0;
                fHistEnergySum[det][strip] -> Fill(sum);
                fHistLeftRight[det][strip] -> Fill(energy1, energy2);
                fHistEnergyPosition[det][strip] -> Fill(pos,sum);
                fHistEnergyPositionAll -> Fill(pos,sum);
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
            //auto cvs2 = MakeCanvas("LeftRight",det);
            auto cvs3 = MakeCanvas("EnergyPosition",det);
            for (auto strip=0; strip<8; ++strip) {
                cvs1 -> cd(strip+1); fHistEnergySum[det][strip] -> Draw();
                //cvs2 -> cd(strip+1); fHistLeftRight[det][strip] -> Draw("colz");
                cvs3 -> cd(strip+1); fHistEnergyPosition[det][strip] -> Draw("colz");
            }
        }
    }
    MakeCanvas("EnergySum",-1,false);
    fHistEnergyPositionAll -> Draw("colz");

    auto fileHist = new TFile(fC2HistFileName,"recreate");
    fHistEnergyPositionAll -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistEnergySum[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistLeftRight[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistEnergyPosition[det][strip] -> Write();

    cout << fC2HistFileName << endl;
}
