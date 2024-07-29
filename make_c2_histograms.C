#include "si_calibration.h"
#include "energy_calibration.h"

void make_c2_histograms(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    /////////////////////////////////////////////////////////////////////
    // 1) Get calibration parameter if calibration is 1 (slope correction
    /////////////////////////////////////////////////////////////////////
    TString calName = "c2_";
    GetC0Parameters();
    GetC1Parameters();
    GetC2Parameters();

    /////////////////////////////////////////////////////////////////////
    // 2) Create histograms
    /////////////////////////////////////////////////////////////////////
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] = MakeHist1(calName+"energy","",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2);
        //fHistEnergy[dss.det][dss.side][dss.strip] = MakeHist1(calName+"energy","",dss.det,dss.side,dss.strip,-1);
    }
    for (auto dss : fStripArrayR)
    {
        fHistEnergySum     [dss.det][dss.side][dss.strip] = MakeHist1(calName+"esum",             "",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2);
        fHistLeftRight     [dss.det][dss.side][dss.strip] = MakeHist2(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1,fNBinE,fBinE1,fBinE2,fNBinE,fBinE1,fBinE2);
        fHistEnergyPosition[dss.det][dss.side][dss.strip] = MakeHist2(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1,fNBinX,fBinX1,fBinX2,fNBinE,fBinE1,1.5*fBinE2);
    }
    fHistEnergyPositionAll = MakeHist2(calName+"rpos",calName+"esum", -1,-1,-1,-1,2*fNBinX,fBinX1,fBinX2,2*fNBinE,fBinE1,1.5*fBinE2);
    fHistEnergyDetector[0] = MakeHist2("det_j",calName+"esum",-1,0,-1,-1,fNumDetectors*fNumJStrips,0,fNumDetectors*fNumJStrips,2*fNBinEFix,fBinE1,1.5*fBinE2);
    fHistEnergyDetector[1] = MakeHist2("det_o",calName+"esum",-1,1,-1,-1,fNumDetectors*fNumOStrips,0,fNumDetectors*fNumOStrips,2*fNBinEFix,fBinE1,1.5*fBinE2);

    /////////////////////////////////////////////////////////////////////
    // 3) Fill histogram from stark event tree
    /////////////////////////////////////////////////////////////////////
    auto fileIn = new TFile(fRecoFileName);
    auto tree = (TTree*) fileIn -> Get("event");
    TClonesArray *array = nullptr;
    tree -> SetBranchAddress("SiChannel",&array);
    auto numEvents = tree -> GetEntries();
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        auto numChannels = array -> GetEntries();
        if (iEvent%20000==0) cout << "Filling raw histogram " << iEvent << " / " << numEvents << " (" << 100*iEvent/numEvents << " %)" << endl;
        if (numChannels!=3) continue;
        for (auto iChannel=0; iChannel<numChannels; ++iChannel)
        {
            auto channel = (LKSiChannel*) array -> At(iChannel);
            auto det = channel -> GetDetID();
            auto side = channel -> GetSide();
            auto strip = channel -> GetStrip();
            auto numStrips = (side==0?fNumJStrips:fNumOStrips);
            if (ContinueRegardingToDataType(det)) continue;
            if (!IsPositionSensitiveStrip(channel))
            {
                double g0 = fC0Parameters[det][side][strip];
                auto energy = channel -> GetEnergy() * g0;
                fHistEnergy[det][side][strip] -> Fill(energy);
                fHistEnergyDetector[side] -> Fill(det*numStrips+strip,energy);
            }
            else
            {
                double g1 = fC1Parameters[det][0][strip][0];
                double g2 = fC1Parameters[det][0][strip][1];
                double b0 = fC2Parameters[det][side][strip][0];
                double b1 = fC2Parameters[det][side][strip][1];
                double b2 = fC2Parameters[det][side][strip][2];
                auto energy1 = channel -> GetEnergy() * g1;
                auto energy2 = channel -> GetEnergy2() * g2;
                if (energy2>0) {
                    auto sum = energy1 + energy2;
                    auto pos = (energy1 - energy2) / sum;
                    sum = sum / (b0 + b1*pos + b2*pos*pos) * f241AmAlphaEnergy1;
                    fHistEnergyDetector[side] -> Fill(det*numStrips+strip,sum);
                    fHistEnergySum     [det][side][strip] -> Fill(sum);
                    fHistLeftRight     [det][side][strip] -> Fill(energy1, energy2);
                    fHistEnergyPosition[det][side][strip] -> Fill(pos,sum);
                    fHistEnergyPositionAll -> Fill(pos,sum);
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 4) Draw examples
    /////////////////////////////////////////////////////////////////////
    for (auto dssGroup : fExampleGroupArrayS)
    {
        int icvs = 1;
        auto cvs0 = dssGroup.MakeGroupCanvas("Energy");
        for (auto dss : dssGroup.array) {
            cvs0 -> cd(icvs++); fHistEnergy[dss.det][dss.side][dss.strip] -> Draw();
        }
        cvs0 -> Modified(); cvs0 -> Update();
    }
    for (auto dssGroup : fExampleGroupArrayR)
    {
        int icvs = 1;
        auto cvs1 = dssGroup.MakeGroupCanvas("EnergySum");
        auto cvs2 = dssGroup.MakeGroupCanvas("LeftRight");
        auto cvs3 = dssGroup.MakeGroupCanvas("EnergyPosition");
        for (auto dss : dssGroup.array)
        {
            cvs1 -> cd(icvs); fHistEnergySum     [dss.det][dss.side][dss.strip] -> Draw();
            cvs2 -> cd(icvs); fHistLeftRight     [dss.det][dss.side][dss.strip] -> Draw("colz");
            cvs3 -> cd(icvs); fHistEnergyPosition[dss.det][dss.side][dss.strip] -> Draw("colz");
            icvs++;
        }
    }
    auto cvs = MakeCanvas("cvs_energy_detector",2);
    cvs -> cd(1); fHistEnergyDetector[0] -> Draw("colz");
    cvs -> cd(2); fHistEnergyDetector[1] -> Draw("colz");
    MakeCanvas("cvs_esum_all",1);
    fHistEnergyPositionAll -> Draw("colz");

    /////////////////////////////////////////////////////////////////////
    // 7) Write histograms
    /////////////////////////////////////////////////////////////////////
    auto fileHist = new TFile(fC2HistFileName,"recreate");
    for (auto dss : fStripArrayS)
    {
        fHistEnergy[dss.det][dss.side][dss.strip] -> Write();
    }
    for (auto dss : fStripArrayR) {
        fHistEnergySum     [dss.det][dss.side][dss.strip] -> Write();
        fHistLeftRight     [dss.det][dss.side][dss.strip] -> Write();
        fHistEnergyPosition[dss.det][dss.side][dss.strip] -> Write();
    }
    fHistEnergyDetector[0] -> Write();
    fHistEnergyDetector[1] -> Write();
    fHistEnergyPositionAll -> Write();
    cout << fC2HistFileName << endl;
}
