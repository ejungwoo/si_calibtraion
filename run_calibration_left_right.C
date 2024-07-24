#include "calibration.h"

void run_calibration_left_right()
{
    TH2D* histLeftRightGate[40][8][2];

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    bool drawFit = true;
    int numGates = 2;
    int chooseGate = 1;

    TString gateFileName = "../data/stark_0253.hist_gated.root";
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    TF1* fitPol = new TF1("fitPol","pol1",0,6000);

    auto fileGate = new TFile(gateFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip)
            for (auto gate=0; gate<numGates; ++gate) {
                histLeftRightGate[det][strip][gate] = (TH2D*) fileGate -> Get(MakeHistName("LeftRight",det,strip,gate));
            }

    //for (auto det=0; det<40; ++det)
    for (auto det : {27})
    {
        TCanvas* cvs = nullptr;
        if (drawFit && det==27) {
            cvs = LKPainter::GetPainter() -> CanvasResize(Form("cvsLeftRightGateFit_%d",det),500,500,0.95);
            cvs -> Divide(4,4,0.01,0.01);
        }
        for (auto strip=0; strip<8; ++strip)
        {
            for (auto gate=0; gate<numGates; ++gate)
            {
                auto hist = histLeftRightGate[det][strip][gate];
                if (gate==chooseGate)
                {
                    auto result = hist -> Fit(fitPol,"Q0N");
                    if (result.Get()!=nullptr)
                        result.Get() -> Print();
                    double intcy = fitPol -> GetParameter(0);
                    double slope = fitPol -> GetParameter(1);
                    double intcx = -intcy/slope;
                    cout << intcx << " " << intcy << endl;
                }
                if (drawFit && det==27) {
                    cvs -> cd(gate+1+strip*2);
                    hist -> Draw("colz");
                }
            }
        }
    }
}
