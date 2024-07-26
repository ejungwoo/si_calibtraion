#include "si_calibration.h"

void run_slope_correction(int run=-1, int chooseGate=1, bool drawExample=true)
{
    MakeRun(run);

    TString foutName = fC1ParFileName;
    ofstream fout(foutName);

    double x0, y0, x1, y1, slope, b, g1, g2;
    TF1* fitPol = new TF1("fitPol","pol1",0,6000);
    TGraph *graph = new TGraph();
    graph -> SetMarkerStyle(20);
    graph -> SetLineColor(kRed);

    auto fileHist = new TFile(fHistFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip) {
            fHistLeftRight[det][strip] = (TH2D*) fileHist -> Get(MakeHistName("LeftRight",det,strip));
            for (auto gate=0; gate<2; ++gate) {
                fHistLeftRightGate[det][strip][gate] = (TH2D*) fileHist -> Get(MakeHistName("LeftRightGated",det,strip,gate));
            }
        }

    for (auto det=0; det<40; ++det)
    {
        cout << det << endl;

        TCanvas* cvs = nullptr;

        bool drawCurrentDet = false;
        if (drawExample && (det==fExampleDet1||det==fExampleDet2))
            drawCurrentDet = true;

        if (drawCurrentDet) {
            cvs = MakeCanvas("LeftRightGateFit",det);
        }

        for (auto strip=0; strip<8; ++strip)
        {
            for (auto gate : {chooseGate})
            {
                auto hist = fHistLeftRight[det][strip];
                auto histGated = fHistLeftRightGate[det][strip][gate];
                if (histGated->GetEntries()<100) {
                    e_warning << histGated -> GetName() << " entries = " << histGated->GetEntries() << endl;
                    g1 = -1;
                    g2 = -1;
                }
                else {
                    histGated -> Fit(fitPol,"Q0N");

                    slope = fitPol -> GetParameter(1);
                    b = fitPol -> GetParameter(0);
                    x0 = 0;
                    y0 = b;
                    x1 = -y0/slope;
                    y1 = 0;

                    g1 = fAlphaEnergy241Am0 / y0;
                    g2 = g1 / slope;
                }

                fout << det << " " << strip << " " << g1 << " " << g2 << endl;

                if (drawCurrentDet)
                {
                    cvs -> cd(strip+1);
                    //histGated -> GetXaxis() -> SetRangeUser(0,x1*1.1);
                    //histGated -> GetYaxis() -> SetRangeUser(0,y0*1.1);
                    //histGated -> Draw("colz");
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
                    lg -> AddEntry((TObject*)nullptr,Form("g1 = %.6f",g1),"");
                    lg -> AddEntry((TObject*)nullptr,Form("g2 = %.6f",g2),"");
                    lg -> Draw();
                }
            }
        }
    }

    cout << foutName << endl;
}
