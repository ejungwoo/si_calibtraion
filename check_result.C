#include "si_calibration.h"
#include "energy_calibration.h"

void run(int det=0)
{
    MakeRun();

    //vector<strip_info> dss_array;

    auto detector = fStark -> GetSiDetector(det);
    auto numJStrips = detector -> GetNumJunctionStrips();
    auto numOStrips = detector -> GetNumOhmicStrips();
    for (auto side=0; side<2; ++side)
    {
        auto numStrips = (side==0?numJStrips:numOStrips);
        for (auto strip=0; strip<numStrips; ++strip) {
            //lk_debug << det << " " << side << " " << strip << endl;
        }
    }

}

TH1D* histEnergy        [fNumDetectors][2][fNumStrips];
TH1D* histEnergySum     [fNumDetectors][2][fNumStrips];
TH2D* histLeftRight     [fNumDetectors][2][fNumStrips];
TH2D* histEnergyPosition[fNumDetectors][2][fNumStrips];
TH2D* histEnergyDetector[2][2];
TH2D* histEnergyDetectorAll[2];

void check_result()
{
    TString calName = "c2_";

    /////////////////////////////////////////////////////////////////////
    MakeRun(199);
    auto file1 = new TFile("data/stark_0199.hist_c2.root","read");
    for (auto dss : fStripArrayS)
    {
        histEnergy[dss.det][dss.side][dss.strip] = (TH1D*) file1 -> Get(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
    }
    for (auto dss : fStripArrayR)
    {
        histEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) file1 -> Get(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
        histLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) file1 -> Get(MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
        histEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) file1 -> Get(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
    }
    histEnergyDetector[0][0] = (TH2D*) file1 -> Get(MakeHistName(calName+"esum","det_j",-1,0,-1,-1));
    histEnergyDetector[0][1] = (TH2D*) file1 -> Get(MakeHistName(calName+"esum","det_o",-1,1,-1,-1));

    /////////////////////////////////////////////////////////////////////
    MakeRun(303);
    auto file2 = new TFile("data/stark_0303.hist_c2.root","read");
    for (auto dss : fStripArrayS)
    {
        histEnergy[dss.det][dss.side][dss.strip] = (TH1D*) file2 -> Get(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
    }
    for (auto dss : fStripArrayR)
    {
        histEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) file2 -> Get(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
        histLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) file2 -> Get(MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
        histEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) file2 -> Get(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
    }
    histEnergyDetector[1][0] = (TH2D*) file2 -> Get(MakeHistName(calName+"esum","det_j",-1,0,-1,-1));
    histEnergyDetector[1][1] = (TH2D*) file2 -> Get(MakeHistName(calName+"esum","det_o",-1,1,-1,-1));
    histEnergyDetectorAll[0] = (TH2D*) histEnergyDetector[0][0] -> Clone("h0");
    histEnergyDetectorAll[1] = (TH2D*) histEnergyDetector[0][1] -> Clone("h1");
    histEnergyDetectorAll[0] -> Add(histEnergyDetector[1][0]);
    histEnergyDetectorAll[1] -> Add(histEnergyDetector[1][1]);

    /////////////////////////////////////////////////////////////////////
    histEnergyDetectorAll[0] -> SetStats(0);
    histEnergyDetectorAll[1] -> SetStats(0);
    histEnergyDetectorAll[0] -> GetXaxis() -> SetTitleOffset(2.7);
    histEnergyDetectorAll[1] -> GetXaxis() -> SetTitleOffset(2.7);
    for (auto det=0; det<fNumDetectors; ++det)
    {
        auto detector = fStark -> GetSiDetector(det);
        auto name = detector -> GetDetTypeName();
        auto ring = detector -> GetLayer();
        TString sring = "dE"; if (ring==1) sring = "E"; if (ring==2) sring = "16E"; 
        histEnergyDetectorAll[0] -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d %s-%s",det,name.Data(),sring.Data()));
        histEnergyDetectorAll[1] -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d %s-%s",det,name.Data(),sring.Data()));
    }
    auto cvs = MakeCanvas("cvs_energy_detector",2);
    auto pad1 = cvs -> cd(1); pad1 -> SetMargin(0.1,0.12,0.18,0.1); histEnergyDetectorAll[0] -> Draw("colz");
    auto pad2 = cvs -> cd(2); pad2 -> SetMargin(0.1,0.12,0.18,0.1); histEnergyDetectorAll[1] -> Draw("colz");

    run(1);
}
