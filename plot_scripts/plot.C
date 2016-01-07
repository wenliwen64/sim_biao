{
    const Float_t PI = 3.1415926;
// 0. Angle-Diff
    TCanvas* can = new TCanvas("can_angle_diff", "can_angle_diff", 1600, 900);
    Float_t small = 1e-12;
    can->Divide(4, 3, small, small);
    TFile* f[12];
    //Read File
    for(Int_t i = 0; i < 12; i++){
	can->cd(i+1);
	gPad->SetLogy();

        TString Name;
	Int_t n_par = (i+1) * 50;
	Name.Form("../sim_biao_root/simulation_1208_%d.root", n_par);
        f[i] = new TFile(Name.Data(), "read");

	// Fuqiang's results
	TH1D* hangle_diff_west_ep = (TH1D*)f[i]->Get("hanglediff_0");
	hangle_diff_west_ep->Scale(1./hangle_diff_west_ep->GetEntries());
	Name.Form("EP orientation diff. with multiplicity %d", n_par);
	hangle_diff_west_ep->SetTitle(Name.Data());
	hangle_diff_west_ep->Draw();

        // Input
	TH1D* hangle_diff_rp = (TH1D*)f[i]->Get("hanglediff_rpv2_rp");
	hangle_diff_rp->Scale(1./hangle_diff_rp->GetEntries());
        hangle_diff_rp->SetLineColor(2);
	hangle_diff_rp->Draw("same");
    }
    can->SaveAs("angle_diff_50to600.eps");
    can->SaveAs("angle_diff_50to600.png");
// 1. MSC 
    TCanvas* can_msc = new TCanvas("can_msc", "can_msc", 1600, 900);
    can_msc->Divide(4, 3, small, small);

    for(Int_t i = 0; i < 12; i++){
        f[i]->Get("msc_os_v2_west_ep");
    }
// 2. Gamma
}
