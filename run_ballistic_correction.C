#include "si_calibration.h"
#include "energy_calibration.h" // for f241AmAlphaEnergy1

void run_ballistic_correction(int run=-1, bool drawExample = true)
{
    MakeRun(run);

    bool usingSymmetricFit = true;

    /////////////////////////////////////////////////////////////////////
    // 1) Get histograms
    /////////////////////////////////////////////////////////////////////
    auto fileHist = new TFile(fC1HistFileName,"read");
    for (auto dss : fStripArrayR) {
        for (auto gate=0; gate<fNumGates; ++gate) {
            fHistEnergyPositionGate[dss.det][dss.side][dss.strip][gate] = (TH2D*) fileHist -> Get(MakeHistName("c1_rpos","c1_esum",dss.det,dss.side,dss.strip,gate));
        }
    }

    /////////////////////////////////////////////////////////////////////
    // 2) Fit 1st order polynomial for left vs right
    /////////////////////////////////////////////////////////////////////
    TF1* fitPol = new TF1("fitPol","pol2",-1,1);
    int det, side, strip;
    double entries, b0, b1, b2;
    auto fpar = new TFile(fC2ParFileName,"recreate");
    auto tree = new TTree("parameters","ballisctic correction parameters b0, b1, b2: x = c1_rpos, c2_esum = c1_esum(x) * source_energy / (b0 + b1*x + b2*x*x)");
    tree -> Branch("det"    ,&det    );
    tree -> Branch("side"   ,&side   );
    tree -> Branch("strip"  ,&strip  );
    tree -> Branch("entries",&entries);
    tree -> Branch("b0"     ,&b0     );
    tree -> Branch("b1"     ,&b1     );
    tree -> Branch("b2"     ,&b2     );
    for (auto dss : fStripArrayR)
    {
        det = dss.det;
        side = dss.side;
        strip = dss.strip;
        auto gate = fChooseGate;
        auto hist = fHistEnergyPositionGate[det][side][strip][gate];
        entries = hist -> GetEntries();
        if (entries<fEntriesCut) {
            e_warning << hist -> GetName() << " entries = " << entries << endl;
            b0 = -1;
            b0 = -1;
            b2 = -1;
        }
        else {
            fitPol -> SetParameter(0,f241AmAlphaEnergy1);
            fitPol -> SetParameter(2,0.2);
            fitPol -> SetParameter(1,0);
            if (usingSymmetricFit)
                fitPol -> FixParameter(1,0);
            hist -> Fit(fitPol,"Q0N");
            b0 = fitPol -> GetParameter(0);
            b1 = fitPol -> GetParameter(1);
            b2 = fitPol -> GetParameter(2);
        }
        //hist -> Draw(); fitPol -> Draw("samel"); cout << b0 << " " << b1 << " " << b2 << endl; return;
        tree -> Fill();
    }
    fpar -> cd();
    tree -> Write();
    fpar -> Close();
    cout << fC2ParFileName << endl;

    /////////////////////////////////////////////////////////////////////
    // 3) Draw examples
    /////////////////////////////////////////////////////////////////////
    GetC2Parameters();
    for (auto dssGroup : fExampleGroupArrayR)
    {
        int icvs = 1;
        auto cvs2 = dssGroup.MakeGroupCanvas("esum_rpos");
        for (auto dss : dssGroup.array)
        {
            cvs2 -> cd(icvs);
            auto det = dss.det;
            auto side = dss.side;
            auto strip = dss.strip;
            auto gate = fChooseGate;
            auto hist = fHistEnergyPositionGate[det][side][strip][gate];
            hist -> Draw("colz");
            TF1* fit = new TF1(Form("fit%d%d%d",det,side,strip),"pol2",-1,1);
            b0 = fC2Parameters[det][side][strip][0];
            b1 = fC2Parameters[det][side][strip][1];
            b2 = fC2Parameters[det][side][strip][2];
            fit -> SetParameter(0,b0);
            fit -> SetParameter(1,b1);
            fit -> SetParameter(2,b2);
            fit -> Draw("samel");
            auto lg = new TLegend(0.45,0.65,0.9,0.88);
            lg -> SetBorderSize(0);
            lg -> SetFillStyle(0);
            lg -> AddEntry((TObject*)nullptr,Form("b0 = %.3f",b0),"");
            lg -> AddEntry((TObject*)nullptr,Form("b1 = %.3f",b1),"");
            lg -> AddEntry((TObject*)nullptr,Form("b2 = %.3f",b2),"");
            //lg -> AddEntry((TObject*)nullptr,Form("y2 = %.1f",fC1Parameters[det][side][strip][7]),"");
            //lg -> AddEntry((TObject*)nullptr,Form("g1 = %.6f",fC1Parameters[det][side][strip][0]),"");
            //lg -> AddEntry((TObject*)nullptr,Form("g2 = %.6f",fC1Parameters[det][side][strip][1]),"");
            lg -> Draw();
            icvs++;
        }
    }
}
