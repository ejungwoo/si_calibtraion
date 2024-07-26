#include "si_calibration.h"

void run_ballistic_correction(int run=-1, int chooseGate = 1, bool drawExample = true)
{
    MakeRun(run);

    bool usingSymmetricFit = true;

    TString foutName = fC2ParFileName;
    ofstream fout(foutName);

    double b0, b1, b2;
    TF1* fitPol;
    if (usingSymmetricFit) fitPol = new TF1("fitPol","[1]*x*x+[0]",-1,1);
    else fitPol = new TF1("fitPol","pol2",-1,1);

    auto fileHist = new TFile(fC1HistFileName,"read");
    for (auto det=0; det<40; ++det)
        for (auto strip=0; strip<8; ++strip) {
            for (auto gate=0; gate<fNumGates; ++gate) {
                fHistEnergyPositionGate[det][strip][gate] = (TH2D*) fileHist -> Get(MakeHistName("C1EnergyPositionGate",det,strip,gate));
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
            cvs = MakeCanvas("C1EnergyPosition",det);
        }

        for (auto strip=0; strip<8; ++strip)
        {
            auto hist = fHistEnergyPositionGate[det][strip][1];

            if (hist->GetEntries()<100) {
                e_warning << hist -> GetName() << " entries = " << hist->GetEntries() << endl;
                b0 = -1;
                b0 = -1;
                b2 = -1;
            }
            else {
                hist -> Fit(fitPol,"Q0N");
                if (usingSymmetricFit) {
                    b0 = fitPol -> GetParameter(0);
                    b1 = 0;
                    b2 = fitPol -> GetParameter(1);
                }
                else {
                    b0 = fitPol -> GetParameter(0);
                    b1 = fitPol -> GetParameter(1);
                    b2 = fitPol -> GetParameter(2);
                }
            }

            fout << det << " " << strip << " " << b0 << " " << b1 << " " << b2 << endl;

            if (drawCurrentDet)
            {
                cvs -> cd(strip+1);
                hist -> Draw("colz");
                fitPol -> DrawClone("samel");
            }
        }
    }

    cout << foutName << endl;
}
