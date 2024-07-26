//////////////////////////////////////////////////////////////////////////////////
      int    fRun    = 253;
const int    fNBinsA = 250;
const double fBinA1  = 0;
const double fBinA2  = 6000;

const int    fNBinsE = 250;
const double fBinE1  = 0;
const double fBinE2  = 10;

const int    fNBinsPos = 200;
const double fBinPos1  = -1;
const double fBinPos2  = 1;

const int    fNumDetectors = 40;
const int    fNumStrips = 8;
const int    fNumGates = 2;

const double fAlphaEnergy241Am0 = 5.486; // 85.2
const double fAlphaEnergy241Am1 = 5.443; // 12.8

const double fExpectedResolution = 0.015;
const int    fExampleDet1 = 14;
const int    fExampleDet2 = 22;
const int    fExampleDet3 = 34;

//////////////////////////////////////////////////////////////////////////////////
TH1D* fHistEnergySum          [fNumDetectors][fNumStrips];
TH2D* fHistLeftRight          [fNumDetectors][fNumStrips];
TH2D* fHistEnergyPosition     [fNumDetectors][fNumStrips];
TH2D* fHistEnergyPositionAll;
TH2D* fHistLeftRightGate      [fNumDetectors][fNumStrips][fNumGates];
TH2D* fHistEnergyPositionGate [fNumDetectors][fNumStrips][fNumGates];
TH2D* fHistLeftRightC1        [fNumDetectors][fNumStrips];
//TH2D* fHistEnergyPositionC1[fNumDetectors][fNumStrips];
TH1D* fHistEnergy             [fNumDetectors][2][fNumStrips];

//////////////////////////////////////////////////////////////////////////////////
TString fRecoFileName;
TString fHistFileName;
//TString fGateFileName;
TString fC1HistFileName;
TString fC2HistFileName;
TString fDummyFileName;
TString fC1ParFileName;
TString fC2ParFileName;

//////////////////////////////////////////////////////////////////////////////////
SKSiArrayPlane* fStark = nullptr;

void MakeRun(int run)
{
    if (run<0) run = fRun;
    fRun = run;
    const char* recoDir = "../data";
    const char* dataDir = "data";
    fRecoFileName   = Form("%s/stark_0%d.reco.root",recoDir,run);
    fHistFileName   = Form("%s/stark_0%d.hist.root",dataDir,run);
    //fGateFileName   = Form("%s/stark_0%d.hist_gated.root",dataDir,run);
    fC1HistFileName = Form("%s/stark_0%d.hist_c1.root",dataDir,run);
    fC1ParFileName  = Form("%s/stark_0%d.c1.txt",dataDir,run);
    fC2HistFileName = Form("%s/stark_0%d.hist_c2.root",dataDir,run);
    fC2ParFileName  = Form("%s/stark_0%d.c2.txt",dataDir,run);
    fDummyFileName  = Form("%s/dummy",dataDir);

    fStark = new SKSiArrayPlane();
    fStark -> AddPar("config_stark.mac");
    fStark -> Init();
}

TString MakeHistName(TString name, int det=-1, int strip=-1, int gate=-1)
{
    if (det<0&&strip<0&&gate<0) return Form("hist%s",name.Data());
    else if (strip<0&&gate<0) return Form("hist%s_%d",name.Data(),det);
    else if (gate<0) return Form("hist%s_%d_%d",name.Data(),det,strip);
    return Form("hist%s_%d_%d_%d",name.Data(),det,strip,gate);
}

TString MakeCanvasName(TString name, int det)
{
    return Form("cvs%s_%d",name.Data(),det);
}

TCanvas *MakeCanvas(TString name, int det=-1, bool divide33=true)
{
    auto cvs = LKPainter::GetPainter() -> CanvasResize(MakeCanvasName(name,det),500,500,0.95);
    if (divide33)
        cvs -> Divide(3,3,0.01,0.01);
    return cvs;
}

bool IsPositionSensitiveStrip(int det, int side=1)
{
    if (det>31) return false;
    if (side==1) return false;
    return true;
}

bool IsPositionSensitiveStrip(LKSiChannel* channel)
{
    return true;
    return IsPositionSensitiveStrip(channel->GetDetID(),channel->GetSide());
}
