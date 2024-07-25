int fNBinsE = 320;
double fBinE1 = 0;
double fBinE2 = 6000;

const int fNumGates = 2;
const double fExpectedResolution = 0.015;

TString fRecoFileName;
TString fHistFileName;
TString fGateFileName;

void SetEnergyBinning(int ne, double e1, double e2)
{
    fNBinsE = ne;
    fBinE1 = e1;
    fBinE2 = e2;
}

void MakeRun(int run)
{
    fRecoFileName = Form("../data/stark_0%d.reco.root",run);
    fHistFileName = Form("../data/stark_0%d.hist.root",run);
    fGateFileName = Form("../data/stark_0%d.hist_gated.root",run);
}

TString MakeHistName(TString name, int det, int strip, int gate=-1)
{
    if (gate<0)
        return Form("hist%s_%d_%d",name.Data(),det,strip);
    return Form("hist%s_%d_%d_%d",name.Data(),det,strip,gate);
}

