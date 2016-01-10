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
	hangle_diff_west_ep->GetYaxis()->CenterTitle();
	hangle_diff_west_ep->GetXaxis()->SetTitle("Angle Diff.");
	hangle_diff_west_ep->GetYaxis()->SetTitle("Normalized prob.");
	gPad->SetTicks(1, 1);
	hangle_diff_west_ep->Draw();

        // Input
	TH1D* hangle_diff_rp = (TH1D*)f[i]->Get("hanglediff_rpv2_rp");
	hangle_diff_rp->Scale(1./hangle_diff_rp->GetEntries());
        hangle_diff_rp->SetLineColor(2);
	hangle_diff_rp->Draw("same");

        TLegend* leg = new TLegend(.1, .7, .3, .9);
	leg->AddEntry(hangle_diff_rp, "input v2=0");
	leg->AddEntry(hangle_diff_west_ep, "reconstructed  v2=0");
	leg->Draw();

    }
    can->SaveAs("angle_diff_50to600.eps");
    can->SaveAs("angle_diff_50to600.png");

// 1. MSC 
    TCanvas* can_msc = new TCanvas("can_msc", "can_msc", 1600, 900);
    can_msc->Divide(4, 3, small, small);

    for(Int_t i = 0; i < 12; i++){
        can_msc->cd(i+1);
	Int_t n_par = (i+1) * 50;
        TProfile* msc_os_rp = (TProfile*)f[i]->Get("msc_os_v2_rp");
        TProfile* msc_ss_rp = (TProfile*)f[i]->Get("msc_ss_v2_rp");

        msc_os_rp->SetMarkerStyle(25); // empty-square
	msc_ss_rp->SetMarkerStyle(25);
	msc_os_rp->SetLineColor(4);
	msc_os_rp->SetMarkerColor(4);
	msc_ss_rp->SetLineColor(2);
	msc_ss_rp->SetMarkerColor(2);

        TString title;
	title.Form("msc input v.s. recon. with n_part_%d", n_par);
	msc_os_rp->SetTitle(title.Data());
	msc_os_rp->GetYaxis()->SetRangeUser(-0.003, 0.003);
	msc_os_rp->GetXaxis()->SetTitle("v2");
	msc_os_rp->GetYaxis()->CenterTitle();
	msc_os_rp->GetYaxis()->SetTitle("msc_{os/ss}");
	gPad->SetTicks(1, 1);
	msc_os_rp->Draw("p");
	msc_ss_rp->Draw("p same");

        TProfile* msc_os_west_ep = (TProfile*)f[i]->Get("msc_os_v2_west_ep");
        TProfile* msc_ss_west_ep = (TProfile*)f[i]->Get("msc_ss_v2_west_ep");

        msc_os_west_ep->SetMarkerStyle(8);
        msc_ss_west_ep->SetMarkerStyle(8);
	msc_os_west_ep->SetLineColor(4);
	msc_os_west_ep->SetMarkerColor(4);
	msc_ss_west_ep->SetLineColor(2);
	msc_ss_west_ep->SetMarkerColor(2);

        msc_os_west_ep->Draw("p same");
        msc_ss_west_ep->Draw("p same");

	TLegend* leg = new TLegend(.4, .14, .65, .32);
	leg->SetBorderSize(0);
	leg->AddEntry(msc_os_rp, "os_real_rp");
	leg->AddEntry(msc_ss_rp, "ss_real_rp");
	leg->AddEntry(msc_os_west_ep, "os_recon_ep");
	leg->AddEntry(msc_ss_west_ep, "ss_recon_ep");
	leg->Draw();
    }
    can_msc->SaveAs("msc_50to60.eps");
    can_msc->SaveAs("msc_50to60.png");

// 2. Gamma
    TCanvas* can_gamma = new TCanvas("can_gamma", "can_gamma", 1600, 900);
    can_gamma->Divide(4, 3, small, small);

    for(Int_t i = 0; i < 12; i++){
        can_gamma->cd(i+1);
	Int_t n_par = (i+1) * 50;
        TProfile* gamma_os_rp = (TProfile*)f[i]->Get("gamma_os_v2_rp");
        TProfile* gamma_ss_rp = (TProfile*)f[i]->Get("gamma_ss_v2_rp");

        gamma_os_rp->SetMarkerStyle(25); // empty-square
	gamma_ss_rp->SetMarkerStyle(25);
	gamma_os_rp->SetLineColor(4);
	gamma_os_rp->SetMarkerColor(4);
	gamma_ss_rp->SetLineColor(2);
	gamma_ss_rp->SetMarkerColor(2);

        TString title;
	title.Form("gamma input v.s. recon. with n_part_%d", n_par);
	gamma_os_rp->SetTitle(title.Data());
	gamma_os_rp->GetYaxis()->SetRangeUser(-0.003, 0.003);
	gamma_os_rp->GetYaxis()->CenterTitle();
	gamma_os_rp->GetXaxis()->SetTitle("v2");
	gamma_os_rp->GetYaxis()->SetTitle("#gamma_{os/ss}");
	gPad->SetTicks(1, 1);
	gamma_os_rp->Draw("p");
	gamma_ss_rp->Draw("p same");

        TProfile* gamma_os_west_ep = (TProfile*)f[i]->Get("gamma_os_v2_west_ep");
        TProfile* gamma_ss_west_ep = (TProfile*)f[i]->Get("gamma_ss_v2_west_ep");

        gamma_os_west_ep->SetMarkerStyle(8);
        gamma_ss_west_ep->SetMarkerStyle(8);
	gamma_os_west_ep->SetLineColor(4);
	gamma_os_west_ep->SetMarkerColor(4);
	gamma_ss_west_ep->SetLineColor(2);
	gamma_ss_west_ep->SetMarkerColor(2);

        gamma_os_west_ep->Draw("p same");
        gamma_ss_west_ep->Draw("p same");

	TLegend* leg = new TLegend(.4, .14, .65, .32);
	leg->SetBorderSize(0);
	leg->AddEntry(gamma_os_rp, "os_real_rp");
	leg->AddEntry(gamma_ss_rp, "ss_real_rp");
	leg->AddEntry(gamma_os_west_ep, "os_recon_ep");
	leg->AddEntry(gamma_ss_west_ep, "ss_recon_ep");
	leg->Draw();
    }
    can_gamma->SaveAs("gamma_50to60.eps");
    can_gamma->SaveAs("gamma_50to60.png");
}
