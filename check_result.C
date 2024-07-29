#include "si_calibration.h"
#include "energy_calibration.h"

TH1D* histEnergy        [fNumDetectors][2][fNumStrips];
TH1D* histEnergySum     [fNumDetectors][2][fNumStrips];
TH2D* histLeftRight     [fNumDetectors][2][fNumStrips];
TH2D* histEnergyPosition[fNumDetectors][2][fNumStrips];
TH2D* histEnergyDetector[2][2];
TH2D* histEnergyDetectorAll[2];

void run(int det=0, bool drawFigures=true)
{
    MakeRun();

    /////////////////////////////////////////////////////////////////////
    // create group
    /////////////////////////////////////////////////////////////////////
    strip_group dssArrayR;
    strip_group dssArrayS;
    auto detector = fStark -> GetSiDetector(det);
    auto numJStrips = detector -> GetNumJunctionStrips();
    auto numOStrips = detector -> GetNumOhmicStrips();
    auto detType = detector -> GetDetType();
    for (auto side=0; side<2; ++side)
    {
        if (IsPositionSensitiveStrip(det,side))
        {
            auto numStrips = (side==0?numJStrips:numOStrips);
            for (auto strip=0; strip<numStrips; ++strip) {
                dssArrayR.array.push_back(strip_info(detType,det,side,strip));
            }
        }
        else {
            auto numStrips = (side==0?numJStrips:numOStrips);
            for (auto strip=0; strip<numStrips; ++strip) {
                dssArrayS.array.push_back(strip_info(detType,det,side,strip));
            }
        }
    }

    /////////////////////////////////////////////////////////////////////
    // create canvas
    /////////////////////////////////////////////////////////////////////
    auto spectrum = new TSpectrum(fNumGates*2);
    auto fitAndDraw = [drawFigures,spectrum](strip_info dss, TH1D* hist)
    {
        auto lg = new TLegend(0.35,0.75,0.80,0.88);
        lg -> SetBorderSize(0);
        lg -> SetFillStyle(0);
        auto numPeaks = spectrum -> Search(hist,5,"goff nodraw");
        double* xPeaks = spectrum -> GetPositionX();
        if (xPeaks[1]<xPeaks[0]) {
            auto xx = xPeaks[1];
            xPeaks[1] = xPeaks[0];
            xPeaks[0] = xx;
        }
        double x1=-999, x2=-999;
        double m1=-999, m2=-999;
        TF1 *fits[2];
        for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
        {
            if (iPeak+1>numPeaks) {
                e_warning << hist -> GetName() << " (" << dss.det << ", " << dss.side << ", " << dss.strip << ") : " << " #peaks = " << numPeaks << " #entries = " << hist -> GetEntries() << endl;
                continue;
            }
            double xPeak = xPeaks[iPeak];
            fits[iPeak] = new TF1(Form("fitGaus_%d_%d_%d_%d",dss.det,dss.side,dss.strip,iPeak),"gaus(0)",0,6000);
            auto fit = fits[iPeak];
            fit -> SetRange(xPeak-5*xPeak*fExpectedResolution,xPeak+5*xPeak*fExpectedResolution);
            fit -> SetParameters(hist->GetBinContent(hist->FindBin(xPeak)),xPeak,xPeak*fExpectedResolution);
            hist -> Fit(fit,"Q0NR");
            auto amp = fit -> GetParameter(0);
            auto mean = fit -> GetParameter(1);
            auto sigma = fit -> GetParameter(2);
            if (iPeak==0) { x1 = mean - 20*sigma; m1 = mean; }
            if (iPeak==1) { x2 = mean + 20*sigma; m2 = mean; }
            fit -> SetRange(mean-2.5*sigma,mean+2.5*sigma);
            lg -> AddEntry((TObject*)nullptr,Form("mean_{%d} = %.3f",iPeak,mean),"");
        }
        if (drawFigures)
        {
            if (x1>-998&&x2>-998)
                hist -> GetXaxis() -> SetRangeUser(x1,x2);
            hist -> Draw();
            for (auto iPeak=0; iPeak<fNumGates; ++iPeak)
                fits[iPeak] -> Draw("samel");
            lg -> Draw("same");
        }
        cout << setw(5) << dss.det << setw(12) << (dss.side==0?"junction":"ohmic") << setw(5) << dss.strip
            << setw(12) << m1 << setw(12) << Form("(%.3f %s)",100*abs(m1-f148GdAlphaEnergy)/f148GdAlphaEnergy,"%")
            << setw(12) << m2 << setw(12) << Form("(%.3f %s)",100*abs(m2-f241AmAlphaEnergy1)/f241AmAlphaEnergy1,"%")
            << endl;
    };

    /////////////////////////////////////////////////////////////////////
    // create canvas
    /////////////////////////////////////////////////////////////////////
    int icvs = 1;
    TCanvas *cvs = nullptr;
    if (drawFigures) {
        cvs = dssArrayS.MakeGroupCanvas("cvsAll",20);
        for (auto dss : dssArrayR.array) { cvs -> cd(icvs++); histEnergyPosition[dss.det][dss.side][dss.strip] -> Draw("colz"); }
    }

    for (auto dss : dssArrayR.array)
    {
        if (drawFigures) cvs -> cd(icvs++);
        auto hist = histEnergySum[dss.det][dss.side][dss.strip];
        fitAndDraw(dss,hist);
    }

    for (auto dss : dssArrayS.array)
    {
        if (drawFigures) cvs -> cd(icvs++);
        auto hist = histEnergy[dss.det][dss.side][dss.strip];
        fitAndDraw(dss,hist);
    }
}

void check_result()
{
    TString calName = "c2_";

    /////////////////////////////////////////////////////////////////////
    // 199
    /////////////////////////////////////////////////////////////////////
    MakeRun(199);
    auto file1 = new TFile(Form("data/stark_0%d.hist_c2.root",fRun),"read");
    for (auto dss : fStripArrayS)
    {
        if (ContinueRegardingToDataType(dss.det)) continue;
        histEnergy[dss.det][dss.side][dss.strip] = (TH1D*) file1 -> Get(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
    }
    for (auto dss : fStripArrayR)
    {
        if (ContinueRegardingToDataType(dss.det)) continue;
        histEnergySum     [dss.det][dss.side][dss.strip] = (TH1D*) file1 -> Get(MakeHistName(calName+"esum",             "",dss.det,dss.side,dss.strip,-1));
        histLeftRight     [dss.det][dss.side][dss.strip] = (TH2D*) file1 -> Get(MakeHistName(calName+"left",calName+"right",dss.det,dss.side,dss.strip,-1));
        histEnergyPosition[dss.det][dss.side][dss.strip] = (TH2D*) file1 -> Get(MakeHistName(calName+"rpos",calName+"esum", dss.det,dss.side,dss.strip,-1));
    }
    histEnergyDetector[0][0] = (TH2D*) file1 -> Get(MakeHistName(calName+"esum","det_j",-1,0,-1,-1));
    histEnergyDetector[0][1] = (TH2D*) file1 -> Get(MakeHistName(calName+"esum","det_o",-1,1,-1,-1));

    /////////////////////////////////////////////////////////////////////
    // 303
    /////////////////////////////////////////////////////////////////////
    MakeRun(303);
    auto file2 = new TFile(Form("data/stark_0%d.hist_c2.root",fRun),"read");
    for (auto dss : fStripArrayS)
    {
        if (ContinueRegardingToDataType(dss.det)) continue;
        histEnergy[dss.det][dss.side][dss.strip] = (TH1D*) file2 -> Get(MakeHistName(calName+"energy","",dss.det,dss.side,dss.strip,-1));
    }
    for (auto dss : fStripArrayR)
    {
        if (ContinueRegardingToDataType(dss.det)) continue;
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
    // merge hit pattern
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
        histEnergyDetectorAll[0] -> GetXaxis() -> SetBinLabel(det*fNumJStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
        histEnergyDetectorAll[1] -> GetXaxis() -> SetBinLabel(det*fNumOStrips+1,Form("%d (%s,%s)",det,name.Data(),sring.Data()));
    }
    auto cvs = MakeCanvas("cvs_energy_detector",2);
    auto pad1 = cvs -> cd(1); pad1 -> SetMargin(0.1,0.12,0.18,0.1); histEnergyDetectorAll[0] -> Draw("colz");
    auto pad2 = cvs -> cd(2); pad2 -> SetMargin(0.1,0.12,0.18,0.1); histEnergyDetectorAll[1] -> Draw("colz");

    run(1);
}
