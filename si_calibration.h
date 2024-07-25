int fRun = 253;

int fNBinsA = 250;
double fBinA1 = 0;
double fBinA2 = 6000;

int fNBinsE = 250;
double fBinE1 = 0;
double fBinE2 = 10;

const int fNumGates = 2;
const double fExpectedResolution = 0.015;

int fExampleDet1 = 14;
int fExampleDet2 = 22;

TH1D* fHistEnergySum[40][8];
TH2D* fHistLeftRight[40][8];
TH2D* fHistEnergyPosition[40][8];
TH2D* fHistLeftRightSlopeCalibrated[40][8];

TString fRecoFileName;
TString fHistFileName;
TString fGateFileName;
TString fSlopeCalibratedHistFileName;

TString fDummyFileName;
TString fSlopeCalibrationParFileName;

void MakeRun(int run)
{
    if (run<0) run = fRun;
    fRun = run;
    const char* recoDir = "../data";
    const char* dataDir = "data";
    fRecoFileName = Form("%s/stark_0%d.reco.root",recoDir,run);
    fHistFileName = Form("%s/stark_0%d.hist.root",dataDir,run);
    fGateFileName = Form("%s/stark_0%d.hist_gated.root",dataDir,run);
    fSlopeCalibratedHistFileName = Form("%s/stark_0%d.hist_slope_calibrated.root",dataDir,run);
    fSlopeCalibrationParFileName = Form("%s/stark_0%d.slope_calibration.txt",dataDir,run);
    fDummyFileName = Form("%s/dummy",dataDir);
}

TString MakeHistName(TString name, int det, int strip, int gate=-1)
{
    if (gate<0)
        return Form("hist%s_%d_%d",name.Data(),det,strip);
    return Form("hist%s_%d_%d_%d",name.Data(),det,strip,gate);
}

TString MakeCanvasName(TString name, int det)
{
    return Form("cvs%s_%d",name.Data(),det);
}

TCanvas *MakeCanvas(TString name, int det)
{
    auto cvs = LKPainter::GetPainter() -> CanvasResize(MakeCanvasName(name,det),500,500,0.95);
    cvs -> Divide(3,3,0.01,0.01);
    return cvs;
}
