#include "si_calibration.h"
#include "energy_calibration.h"

void run_slope_correction(int run=-1, bool drawExample=true)
{
    MakeRun(run);

    /////////////////////////////////////////////////////////////////////
    // 1) Get histograms
    /////////////////////////////////////////////////////////////////////
    auto fileHist = new TFile(fHistFileName,"read");
    for (auto dss : fStripArrayR)
    {
        fHistLeftRight[dss.det][dss.side][dss.strip] = (TH2D*) fileHist -> Get(MakeHistName("left","right",dss.det,dss.side,dss.strip));
        for (auto gate=0; gate<fNumGates; ++gate) {
            fHistLeftRightGate[dss.det][dss.side][dss.strip][gate] = (TH2D*) fileHist -> Get(MakeHistName("left","right",dss.det,dss.side,dss.strip,gate));
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 2) Fit 1st order polynomial for left vs right
    /////////////////////////////////////////////////////////////////////
    int det, side, strip;
    double entries, slope, b, x1, y1, x2, y2, g1, g2;
    auto fout = new TFile(fC1ParFileName,"recreate");
    auto tree = new TTree("parameters","slope correction parameters g1, g2: energy_left  = adc_left * g1, energy_right = adc_right * g2");
    tree -> Branch("det"    ,&det    );
    tree -> Branch("side"   ,&side   );
    tree -> Branch("strip"  ,&strip  );
    tree -> Branch("entries",&entries);
    tree -> Branch("x1"     ,&x1     );
    tree -> Branch("y1"     ,&y1     );
    tree -> Branch("x2"     ,&x2     );
    tree -> Branch("y2"     ,&y2     );
    tree -> Branch("b"      ,&b      );
    tree -> Branch("slope"  ,&slope  );
    tree -> Branch("g1"     ,&g1     );
    tree -> Branch("g2"     ,&g2     );
    TF1* fitPol = new TF1("fitPol","pol1",0,6000);
    for (auto dss : fStripArrayR)
    {
        det = dss.det;
        side = dss.side;
        strip = dss.strip;
        auto gate = fChooseGate;
        auto hist = fHistLeftRight[det][side][strip];
        auto histGated = fHistLeftRightGate[det][side][strip][gate];
        entries = histGated->GetEntries();
        if (entries<fEntriesCut) {
            e_warning << histGated -> GetName() << " entries = " << entries << endl;
            slope = 0;
            b = 0;
            x1 = 0;
            y1 = 0;
            x2 = 0;
            y2 = 0;
            g1=-1;
            g2=-1;
        }
        else {
            histGated -> Fit(fitPol,"Q0N");
            slope = fitPol -> GetParameter(1);
            b = fitPol -> GetParameter(0);
            x1 = 0;
            y1 = b;
            x2 = -y1/slope;
            y2 = 0;
            g1 = f241AmAlphaEnergy1 / y1;
            g2 = -g1 / slope;
        }
        tree -> Fill();
    }
    fout -> cd();
    tree -> Write();
    fout -> Close();
    cout << fC1ParFileName << endl;

    /////////////////////////////////////////////////////////////////////
    // 3) Draw examples
    /////////////////////////////////////////////////////////////////////
    GetC1Parameters();
    auto graph = new TGraph();
    graph -> SetMarkerStyle(20);
    graph -> SetLineColor(kRed);
    for (auto dssGroup : fExampleGroupArrayR)
    {
        int icvs = 1;
        auto cvs2 = dssGroup.MakeGroupCanvas("LeftRight");
        for (auto dss : dssGroup.array)
        {
            auto det = dss.det;
            auto side = dss.side;
            auto strip = dss.strip;
            cvs2 -> cd(icvs);
            fHistLeftRight[dss.det][dss.side][dss.strip] -> Draw("colz");
            graph -> SetPoint(0,fC1Parameters[det][side][strip][4], fC1Parameters[det][side][strip][5]);
            graph -> SetPoint(1,fC1Parameters[det][side][strip][6], fC1Parameters[det][side][strip][7]);
            graph -> DrawClone("samelp");
            auto lg = new TLegend(0.4,0.55,0.9,0.88);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            lg -> AddEntry((TObject*)nullptr,Form("slope = %.3f",fC1Parameters[det][side][strip][3]),"");
            lg -> AddEntry((TObject*)nullptr,Form("x1 = %.1f",fC1Parameters[det][side][strip][4]),"");
            lg -> AddEntry((TObject*)nullptr,Form("y1 = %.1f",fC1Parameters[det][side][strip][5]),"");
            lg -> AddEntry((TObject*)nullptr,Form("x1 = %.1f",fC1Parameters[det][side][strip][6]),"");
            lg -> AddEntry((TObject*)nullptr,Form("y2 = %.1f",fC1Parameters[det][side][strip][7]),"");
            lg -> AddEntry((TObject*)nullptr,Form("g1 = %.6f",fC1Parameters[det][side][strip][0]),"");
            lg -> AddEntry((TObject*)nullptr,Form("g2 = %.6f",fC1Parameters[det][side][strip][1]),"");
            lg -> Draw();
            icvs++;
        }
    }
}
