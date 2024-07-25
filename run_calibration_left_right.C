#include "calibration.h"
#include "LKLogger.h"

void run_calibration_left_right()
{
    MakeRun(253);

    TH2D* histLeftRightGate[40][8][2];
    TH2D* histLeftRight[40][8];

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    bool drawFit = true;
    int numGates = 2;
    int chooseGate = 1;

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    double x0, y0, x1, y1, slope;
    TF1* fitPol = new TF1("fitPol","pol1",0,6000);
    TGraph *graph = new TGraph();
    graph -> SetMarkerStyle(20);
    //graph -> SetMarkerSize(0.8);
    graph -> SetLineColor(kRed);

    auto fileHist = new TFile(fHistFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip)
            histLeftRight[det][strip] = (TH2D*) fileHist -> Get(MakeHistName("LeftRight",det,strip));

    auto fileGate = new TFile(fGateFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip)
            for (auto gate=0; gate<numGates; ++gate) {
                histLeftRightGate[det][strip][gate] = (TH2D*) fileGate -> Get(MakeHistName("LeftRightGated",det,strip,gate));
            }

    //for (auto det=0; det<40; ++det)
    for (auto det : {27})
    {
        TCanvas* cvs = nullptr;
        if (drawFit && det==27) {
            cvs = LKPainter::GetPainter() -> CanvasResize(Form("cvsLeftRightGateFit_%d",det),500,500,0.95);
            cvs -> Divide(3,3,0.01,0.01);
        }
        for (auto strip=0; strip<8; ++strip)
        {
            for (auto gate=0; gate<numGates; ++gate)
            {
                auto hist = histLeftRight[det][strip];
                auto histGated = histLeftRightGate[det][strip][gate];
                if (gate==chooseGate)
                {
                    if (histGated->GetEntries()<100) {
                        e_error << histGated -> GetName() << " entries = " << histGated->GetEntries() << endl;
                        continue;
                    }

                    histGated -> Fit(fitPol,"Q0N");
                    slope = fitPol -> GetParameter(1);
                    x0 = 0;
                    y0 = fitPol -> GetParameter(0);
                    x1 = -y0/slope;
                    y1 = 0;

                    if (drawFit) {
                        cvs -> cd(strip+1);
                        hist -> GetXaxis() -> SetRangeUser(0,x1*1.1);
                        hist -> GetYaxis() -> SetRangeUser(0,y0*1.1);
                        hist -> Draw("colz");
                        graph -> SetPoint(0,x0,y0);
                        graph -> SetPoint(1,x1,y1);
                        graph -> DrawClone("samepl");
                        auto lg = new TLegend(0.4,0.55,0.9,0.88);
                        lg -> SetBorderSize(0);
                        lg -> SetFillStyle(0);
                        lg -> AddEntry((TObject*)nullptr,Form("slope = %.3f",slope),"");
                        lg -> AddEntry((TObject*)nullptr,Form("x0 = %.1f",x0),"");
                        lg -> AddEntry((TObject*)nullptr,Form("y0 = %.1f",y0),"");
                        lg -> AddEntry((TObject*)nullptr,Form("x1 = %.1f",x1),"");
                        lg -> AddEntry((TObject*)nullptr,Form("y1 = %.1f",y1),"");
                        lg -> Draw();
                    }
                }
            }
        }
    }
}
