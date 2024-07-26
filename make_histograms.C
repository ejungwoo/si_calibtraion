#include "si_calibration.h"

void make_histograms(int calibration=1, int run=-1, bool drawExample=true)
{
    MakeRun(run);

    int ne = fNBinsA;
    double e1 = fBinA1;
    double e2 = fBinA2;

    double gArray[40][8][2];
    TString calibrationName;
    if (calibration==1)
    {
        calibrationName = "C1";
        ifstream fcalibration(fC1ParFileName);
        int det, strip;
        double g1, g2;
        while (fcalibration >> det >> strip >> g1 >> g2)
        {
            gArray[det][strip][0] = abs(g1);
            gArray[det][strip][1] = abs(g2);
        }
        fcalibration.close();
        ne = fNBinsE;
        e1 = fBinE1;
        e2 = fBinE2;
    }

    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);

    //auto numDetectors = fStark -> GetNumSiDetectors();
    //for (auto det=0; det<numDetectors; ++det) {
    //    auto detector = fStark -> GetSiDetector(det);
        //lk_debug << detector -> GetName() << endl;
    //}

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
            for (auto gate=0; gate<fNumGates; ++gate) {
                fHistLeftRightGate[det][strip][gate] = new TH2D(MakeHistName(calibrationName+"LeftRightGated",det,strip,gate),title+";left;right",ne,e1,e2,ne,e1,e2);
                fHistLeftRightGate[det][strip][gate] -> SetStats(0);
                fHistEnergyPositionGate[det][strip][gate] = new TH2D(MakeHistName(calibrationName+"EnergyPositionGate",det,strip,gate),title+";position;energy sum",200,-1,1,ne,e1,1.5*e2);
                fHistEnergyPositionGate[det][strip][gate] -> SetStats(0);
            }
            for (auto side=0; side<2; ++side) {
                fHistEnergy[det][side][strip] = new TH1D(MakeHistName(calibrationName+"Energy",det,side,strip),title+";energy;",ne,e1,e2);
                fHistEnergy[det][side][strip] -> SetStats(0);
            }
        }
    }

    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0)
            cout << "Filling raw histogram " << iEvent << " / " << numEvents << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto strip = channel -> GetStrip();
            if (!IsPositionSensitiveStrip(channel))
            {
                double c1 = 1.;
                auto side = channel -> GetSide();
                auto energy = channel -> GetEnergy() * c1;
                fHistEnergy[det][side][strip] -> Fill(energy);
            }
            else {
                double g1 = 1.;
                double g2 = 1.;
                if (calibration==1) {
                    g1 = gArray[det][strip][0];
                    g2 = gArray[det][strip][1];
                }
                auto energy1 = channel -> GetEnergy() * g1;
                auto energy2 = channel -> GetEnergy2() * g2;
                if (energy2>0) {
                    auto sum = energy1 + energy2;
                    auto pos = (energy1 - energy2) / sum;
                    fHistEnergySum[det][strip] -> Fill(sum);
                    fHistLeftRight[det][strip] -> Fill(energy1, energy2);
                    fHistEnergyPosition[det][strip] -> Fill(pos,sum);
                }
            }
        }
    }

    for (auto det=0; det<40; ++det)
    {
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

    double gatingRange[40][8][fNumGates][2] = {0}; // det, strip, gate, range
    auto spectrum = new TSpectrum(fNumGates*2);
    TF1 *fitGaus = new TF1("fitGaus","gaus(0)",0,6000);

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
            auto hist = fHistEnergySum[det][strip];
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

    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0)
            cout << "Filling gated histogram " << iEvent << " / " << numEvents << endl;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto strip = channel -> GetStrip();
            if (!IsPositionSensitiveStrip(channel))
            {
            }
            else {
                double g1 = 1.;
                double g2 = 1.;
                if (calibration==1) {
                    g1 = gArray[det][strip][0];
                    g2 = gArray[det][strip][1];
                }
                auto energy1 = channel -> GetEnergy() * g1;
                auto energy2 = channel -> GetEnergy2() * g2;
                if (energy2>0) {
                    auto sum = energy1 + energy2;
                    auto pos = (energy1 - energy2) / sum;
                    for (auto gate=0; gate<fNumGates; ++gate) {
                        if (sum>gatingRange[det][strip][gate][0] && sum<gatingRange[det][strip][gate][1]) {
                            fHistLeftRightGate[det][strip][gate] -> Fill(energy1, energy2);
                            fHistEnergyPositionGate[det][strip][gate] -> Fill(pos,sum);
                        }
                    }
                }
            }
        }
    }

    TString foutName = fHistFileName;
    if (calibration==1)
        foutName = fC1HistFileName;

    auto fileHist = new TFile(foutName,"recreate");
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistEnergySum[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistLeftRight[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) fHistEnergyPosition[det][strip] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) for (auto gate=0; gate<fNumGates; ++gate) fHistLeftRightGate[det][strip][gate] -> Write();
    for (auto det=0; det<40; ++det) for (auto strip=0; strip<8; ++strip) for (auto gate=0; gate<fNumGates; ++gate) fHistEnergyPositionGate[det][strip][gate] -> Write();
    for (auto det=0; det<40; ++det) for (auto side=0; side<2; ++side) for (auto strip=0; strip<8; ++strip) fHistEnergy[det][side][strip] -> Write();

    cout << foutName << endl;
}
