#ifdef __CLING__
#pragma cling optimize(0)
#endif
void ptDist()
{
//=========Macro generated from canvas: cPN/
//=========  (Sat Aug 26 16:27:19 2023) by ROOT version 6.28/04
   TCanvas *cPN = new TCanvas("cPN", "",78,53,900,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cPN->SetHighLightColor(2);
   cPN->Range(-0.3333333,0.05424478,2.047619,5.426955);
   cPN->SetFillColor(0);
   cPN->SetBorderMode(0);
   cPN->SetBorderSize(2);
   cPN->SetLogy();
   cPN->SetLeftMargin(0.14);
   cPN->SetRightMargin(0.02);
   cPN->SetTopMargin(0.06);
   cPN->SetBottomMargin(0.12);
   cPN->SetFrameBorderMode(0);
   cPN->SetFrameBorderMode(0);
   
   TH1F *hPt_all__13 = new TH1F("hPt_all__13","ALICE Run-4 Simulation: Pb#minusPb UPC at #sqrt{#it{s}_{NN}} = 5.02 TeV, J/#psi and #psi' #rightarrow e^{+}e^{-}",40,0,2);
   hPt_all__13->SetBinContent(1,11775.76);
   hPt_all__13->SetBinContent(2,28375.39);
   hPt_all__13->SetBinContent(3,28850.54);
   hPt_all__13->SetBinContent(4,23025.28);
   hPt_all__13->SetBinContent(5,17675.99);
   hPt_all__13->SetBinContent(6,14314.28);
   hPt_all__13->SetBinContent(7,11396.93);
   hPt_all__13->SetBinContent(8,9384.872);
   hPt_all__13->SetBinContent(9,8332.949);
   hPt_all__13->SetBinContent(10,7318.751);
   hPt_all__13->SetBinContent(11,7004.7);
   hPt_all__13->SetBinContent(12,5655.693);
   hPt_all__13->SetBinContent(13,4913.107);
   hPt_all__13->SetBinContent(14,4640.273);
   hPt_all__13->SetBinContent(15,3665.831);
   hPt_all__13->SetBinContent(16,2862.284);
   hPt_all__13->SetBinContent(17,2375.11);
   hPt_all__13->SetBinContent(18,1858.358);
   hPt_all__13->SetBinContent(19,1284.538);
   hPt_all__13->SetBinContent(20,1052.245);
   hPt_all__13->SetBinContent(21,831.5363);
   hPt_all__13->SetBinContent(22,534.9962);
   hPt_all__13->SetBinContent(23,365.6916);
   hPt_all__13->SetBinContent(24,274.4707);
   hPt_all__13->SetBinContent(25,211.4299);
   hPt_all__13->SetBinContent(26,118.0858);
   hPt_all__13->SetBinContent(27,94.52309);
   hPt_all__13->SetBinContent(28,34.4753);
   hPt_all__13->SetBinContent(29,56.60955);
   hPt_all__13->SetBinContent(30,28.56553);
   hPt_all__13->SetBinContent(31,28.04402);
   hPt_all__13->SetBinContent(32,9.348006);
   hPt_all__13->SetBinContent(33,9.348006);
   hPt_all__13->SetBinError(1,305.0658);
   hPt_all__13->SetBinError(2,474.5832);
   hPt_all__13->SetBinError(3,478.4046);
   hPt_all__13->SetBinError(4,428.367);
   hPt_all__13->SetBinError(5,376.7174);
   hPt_all__13->SetBinError(6,340.5176);
   hPt_all__13->SetBinError(7,305.5931);
   hPt_all__13->SetBinError(8,281.2009);
   hPt_all__13->SetBinError(9,265.5218);
   hPt_all__13->SetBinError(10,251.0523);
   hPt_all__13->SetBinError(11,247.2018);
   hPt_all__13->SetBinError(12,223.7298);
   hPt_all__13->SetBinError(13,208.1149);
   hPt_all__13->SetBinError(14,203.2043);
   hPt_all__13->SetBinError(15,180.596);
   hPt_all__13->SetBinError(16,159.8655);
   hPt_all__13->SetBinError(17,145.8295);
   hPt_all__13->SetBinError(18,128.7584);
   hPt_all__13->SetBinError(19,106.7153);
   hPt_all__13->SetBinError(20,97.60866);
   hPt_all__13->SetBinError(21,85.50859);
   hPt_all__13->SetBinError(22,69.57281);
   hPt_all__13->SetBinError(23,57.50677);
   hPt_all__13->SetBinError(24,49.55302);
   hPt_all__13->SetBinError(25,43.66254);
   hPt_all__13->SetBinError(26,32.75442);
   hPt_all__13->SetBinError(27,29.57019);
   hPt_all__13->SetBinError(28,16.93094);
   hPt_all__13->SetBinError(29,22.90378);
   hPt_all__13->SetBinError(30,16.19962);
   hPt_all__13->SetBinError(31,16.19122);
   hPt_all__13->SetBinError(32,9.348006);
   hPt_all__13->SetBinError(33,9.348006);
   hPt_all__13->SetMinimum(5);
   hPt_all__13->SetMaximum(127230.9);
   hPt_all__13->SetEntries(26141);
   hPt_all__13->SetStats(0);
   hPt_all__13->SetLineWidth(2);
   hPt_all__13->SetMarkerStyle(24);
   hPt_all__13->SetMarkerSize(1.2);
   hPt_all__13->GetXaxis()->SetTitle("#it{p}_{T,cl pair} (GeV/#it{c})");
   hPt_all__13->GetXaxis()->SetLabelFont(42);
   hPt_all__13->GetXaxis()->SetLabelSize(0.05);
   hPt_all__13->GetXaxis()->SetTitleSize(0.05);
   hPt_all__13->GetXaxis()->SetTitleOffset(1.05);
   hPt_all__13->GetXaxis()->SetTitleFont(42);
   hPt_all__13->GetYaxis()->SetTitle("counts per 50 MeV/#it{c}");
   hPt_all__13->GetYaxis()->SetNdivisions(3000510);
   hPt_all__13->GetYaxis()->SetLabelFont(42);
   hPt_all__13->GetYaxis()->SetLabelSize(0.05);
   hPt_all__13->GetYaxis()->SetTitleSize(0.05);
   hPt_all__13->GetYaxis()->SetTitleOffset(1.25);
   hPt_all__13->GetYaxis()->SetTitleFont(42);
   hPt_all__13->GetZaxis()->SetLabelFont(42);
   hPt_all__13->GetZaxis()->SetTitleOffset(1);
   hPt_all__13->GetZaxis()->SetTitleFont(42);
   hPt_all__13->Draw("E");
   
   TH1F *hPt_cohJpsi__14 = new TH1F("hPt_cohJpsi__14","coh J/#psi #rightarrow e^{+}e^{-}",40,0,2);
   hPt_cohJpsi__14->SetBinContent(1,10509.63);
   hPt_cohJpsi__14->SetBinContent(2,25196.16);
   hPt_cohJpsi__14->SetBinContent(3,23571.37);
   hPt_cohJpsi__14->SetBinContent(4,16747.25);
   hPt_cohJpsi__14->SetBinContent(5,10248.08);
   hPt_cohJpsi__14->SetBinContent(6,6277.245);
   hPt_cohJpsi__14->SetBinContent(7,3740.984);
   hPt_cohJpsi__14->SetBinContent(8,1933.899);
   hPt_cohJpsi__14->SetBinContent(9,966.9493);
   hPt_cohJpsi__14->SetBinContent(10,554.807);
   hPt_cohJpsi__14->SetBinContent(11,142.6647);
   hPt_cohJpsi__14->SetBinContent(12,63.40651);
   hPt_cohJpsi__14->SetBinContent(13,39.62907);
   hPt_cohJpsi__14->SetBinContent(15,7.925814);
   hPt_cohJpsi__14->SetBinError(1,288.6128);
   hPt_cohJpsi__14->SetBinError(2,446.8782);
   hPt_cohJpsi__14->SetBinError(3,432.2295);
   hPt_cohJpsi__14->SetBinError(4,364.3289);
   hPt_cohJpsi__14->SetBinError(5,284.9989);
   hPt_cohJpsi__14->SetBinError(6,223.0522);
   hPt_cohJpsi__14->SetBinError(7,172.1928);
   hPt_cohJpsi__14->SetBinError(8,123.8052);
   hPt_cohJpsi__14->SetBinError(9,87.54348);
   hPt_cohJpsi__14->SetBinError(10,66.31212);
   hPt_cohJpsi__14->SetBinError(11,33.62638);
   hPt_cohJpsi__14->SetBinError(12,22.41759);
   hPt_cohJpsi__14->SetBinError(13,17.72266);
   hPt_cohJpsi__14->SetBinError(15,7.925814);
   hPt_cohJpsi__14->SetMinimum(1);
   hPt_cohJpsi__14->SetMaximum(27778.77);
   hPt_cohJpsi__14->SetEntries(12617);
   hPt_cohJpsi__14->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   hPt_cohJpsi__14->SetLineColor(ci);
   hPt_cohJpsi__14->SetLineWidth(3);
   hPt_cohJpsi__14->GetXaxis()->SetTitle("#it{m}_{supcl pair} [GeV/#it{c}^{2}]");
   hPt_cohJpsi__14->GetXaxis()->SetDecimals();
   hPt_cohJpsi__14->GetXaxis()->SetLabelFont(42);
   hPt_cohJpsi__14->GetXaxis()->SetTitleOffset(1.2);
   hPt_cohJpsi__14->GetXaxis()->SetTitleFont(42);
   hPt_cohJpsi__14->GetYaxis()->SetTitle("Counts per 80 MeV [-]");
   hPt_cohJpsi__14->GetYaxis()->SetDecimals();
   hPt_cohJpsi__14->GetYaxis()->SetNdivisions(3000510);
   hPt_cohJpsi__14->GetYaxis()->SetLabelFont(42);
   hPt_cohJpsi__14->GetYaxis()->SetTitleFont(42);
   hPt_cohJpsi__14->GetZaxis()->SetLabelFont(42);
   hPt_cohJpsi__14->GetZaxis()->SetTitleOffset(1);
   hPt_cohJpsi__14->GetZaxis()->SetTitleFont(42);
   hPt_cohJpsi__14->Draw("HIST SAME");
   
   TH1F *hPt_incJpsi__15 = new TH1F("hPt_incJpsi__15","inc J/#psi #rightarrow e^{+}e^{-}",40,0,2);
   hPt_incJpsi__15->SetBinContent(1,738.4925);
   hPt_incJpsi__15->SetBinContent(2,2112.649);
   hPt_incJpsi__15->SetBinContent(3,3421.37);
   hPt_incJpsi__15->SetBinContent(4,4262.691);
   hPt_incJpsi__15->SetBinContent(5,5188.144);
   hPt_incJpsi__15->SetBinContent(6,5823.808);
   hPt_incJpsi__15->SetBinContent(7,5730.328);
   hPt_incJpsi__15->SetBinContent(8,6048.16);
   hPt_incJpsi__15->SetBinContent(9,5982.724);
   hPt_incJpsi__15->SetBinContent(10,5749.024);
   hPt_incJpsi__15->SetBinContent(11,5954.68);
   hPt_incJpsi__15->SetBinContent(12,5029.228);
   hPt_incJpsi__15->SetBinContent(13,4346.823);
   hPt_incJpsi__15->SetBinContent(14,4215.951);
   hPt_incJpsi__15->SetBinContent(15,3327.89);
   hPt_incJpsi__15->SetBinContent(16,2626.79);
   hPt_incJpsi__15->SetBinContent(17,2187.433);
   hPt_incJpsi__15->SetBinContent(18,1701.337);
   hPt_incJpsi__15->SetBinContent(19,1159.153);
   hPt_incJpsi__15->SetBinContent(20,990.8887);
   hPt_incJpsi__15->SetBinContent(21,738.4925);
   hPt_incJpsi__15->SetBinContent(22,504.7923);
   hPt_incJpsi__15->SetBinContent(23,345.8762);
   hPt_incJpsi__15->SetBinContent(24,252.3962);
   hPt_incJpsi__15->SetBinContent(25,196.3081);
   hPt_incJpsi__15->SetBinContent(26,112.1761);
   hPt_incJpsi__15->SetBinContent(27,93.48006);
   hPt_incJpsi__15->SetBinContent(28,28.04402);
   hPt_incJpsi__15->SetBinContent(29,56.08804);
   hPt_incJpsi__15->SetBinContent(30,28.04402);
   hPt_incJpsi__15->SetBinContent(31,28.04402);
   hPt_incJpsi__15->SetBinContent(32,9.348006);
   hPt_incJpsi__15->SetBinContent(33,9.348006);
   hPt_incJpsi__15->SetBinError(1,83.0869);
   hPt_incJpsi__15->SetBinError(2,140.5313);
   hPt_incJpsi__15->SetBinError(3,178.8379);
   hPt_incJpsi__15->SetBinError(4,199.6188);
   hPt_incJpsi__15->SetBinError(5,220.2244);
   hPt_incJpsi__15->SetBinError(6,233.3259);
   hPt_incJpsi__15->SetBinError(7,231.4458);
   hPt_incJpsi__15->SetBinError(8,237.7777);
   hPt_incJpsi__15->SetBinError(9,236.4879);
   hPt_incJpsi__15->SetBinError(10,231.823);
   hPt_incJpsi__15->SetBinError(11,235.933);
   hPt_incJpsi__15->SetBinError(12,216.8254);
   hPt_incJpsi__15->SetBinError(13,201.5791);
   hPt_incJpsi__15->SetBinError(14,198.5214);
   hPt_incJpsi__15->SetBinError(15,176.3778);
   hPt_incJpsi__15->SetBinError(16,156.7011);
   hPt_incJpsi__15->SetBinError(17,142.997);
   hPt_incJpsi__15->SetBinError(18,126.1115);
   hPt_incJpsi__15->SetBinError(19,104.095);
   hPt_incJpsi__15->SetBinError(20,96.24361);
   hPt_incJpsi__15->SetBinError(21,83.0869);
   hPt_incJpsi__15->SetBinError(22,68.69354);
   hPt_incJpsi__15->SetBinError(23,56.8617);
   hPt_incJpsi__15->SetBinError(24,48.57366);
   hPt_incJpsi__15->SetBinError(25,42.83795);
   hPt_incJpsi__15->SetBinError(26,32.38244);
   hPt_incJpsi__15->SetBinError(27,29.56099);
   hPt_incJpsi__15->SetBinError(28,16.19122);
   hPt_incJpsi__15->SetBinError(29,22.89785);
   hPt_incJpsi__15->SetBinError(30,16.19122);
   hPt_incJpsi__15->SetBinError(31,16.19122);
   hPt_incJpsi__15->SetBinError(32,9.348006);
   hPt_incJpsi__15->SetBinError(33,9.348006);
   hPt_incJpsi__15->SetMinimum(1);
   hPt_incJpsi__15->SetMaximum(6668.097);
   hPt_incJpsi__15->SetEntries(8451);
   hPt_incJpsi__15->SetStats(0);

   ci = TColor::GetColor("#ff0000");
   hPt_incJpsi__15->SetLineColor(ci);
   hPt_incJpsi__15->SetLineWidth(3);
   hPt_incJpsi__15->GetXaxis()->SetTitle("#it{m}_{supcl pair} [GeV/#it{c}^{2}]");
   hPt_incJpsi__15->GetXaxis()->SetDecimals();
   hPt_incJpsi__15->GetXaxis()->SetLabelFont(42);
   hPt_incJpsi__15->GetXaxis()->SetTitleOffset(1.2);
   hPt_incJpsi__15->GetXaxis()->SetTitleFont(42);
   hPt_incJpsi__15->GetYaxis()->SetTitle("Counts per 80 MeV [-]");
   hPt_incJpsi__15->GetYaxis()->SetDecimals();
   hPt_incJpsi__15->GetYaxis()->SetNdivisions(3000510);
   hPt_incJpsi__15->GetYaxis()->SetLabelFont(42);
   hPt_incJpsi__15->GetYaxis()->SetTitleFont(42);
   hPt_incJpsi__15->GetZaxis()->SetLabelFont(42);
   hPt_incJpsi__15->GetZaxis()->SetTitleOffset(1);
   hPt_incJpsi__15->GetZaxis()->SetTitleFont(42);
   hPt_incJpsi__15->Draw("HIST SAME");
   
   TH1F *hPt_cohFD__16 = new TH1F("hPt_cohFD__16","coh #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",40,0,2);
   hPt_cohFD__16->SetBinContent(1,394.2029);
   hPt_cohFD__16->SetBinContent(2,828.9855);
   hPt_cohFD__16->SetBinContent(3,1362.319);
   hPt_cohFD__16->SetBinContent(4,1431.884);
   hPt_cohFD__16->SetBinContent(5,1600);
   hPt_cohFD__16->SetBinContent(6,1379.71);
   hPt_cohFD__16->SetBinContent(7,1171.015);
   hPt_cohFD__16->SetBinContent(8,701.4493);
   hPt_cohFD__16->SetBinContent(9,550.7246);
   hPt_cohFD__16->SetBinContent(10,284.058);
   hPt_cohFD__16->SetBinContent(11,173.913);
   hPt_cohFD__16->SetBinContent(12,52.17391);
   hPt_cohFD__16->SetBinContent(13,46.37681);
   hPt_cohFD__16->SetBinContent(14,17.3913);
   hPt_cohFD__16->SetBinContent(16,5.797101);
   hPt_cohFD__16->SetBinError(1,47.80412);
   hPt_cohFD__16->SetBinError(2,69.32325);
   hPt_cohFD__16->SetBinError(3,88.86788);
   hPt_cohFD__16->SetBinError(4,91.1086);
   hPt_cohFD__16->SetBinError(5,96.30868);
   hPt_cohFD__16->SetBinError(6,89.43333);
   hPt_cohFD__16->SetBinError(7,82.39229);
   hPt_cohFD__16->SetBinError(8,63.76812);
   hPt_cohFD__16->SetBinError(9,56.50316);
   hPt_cohFD__16->SetBinError(10,40.57971);
   hPt_cohFD__16->SetBinError(11,31.75203);
   hPt_cohFD__16->SetBinError(12,17.3913);
   hPt_cohFD__16->SetBinError(13,16.39668);
   hPt_cohFD__16->SetBinError(14,10.04087);
   hPt_cohFD__16->SetBinError(16,5.797101);
   hPt_cohFD__16->SetMinimum(1);
   hPt_cohFD__16->SetMaximum(1764);
   hPt_cohFD__16->SetEntries(1725);
   hPt_cohFD__16->SetStats(0);

   ci = TColor::GetColor("#0000ff");
   hPt_cohFD__16->SetLineColor(ci);
   hPt_cohFD__16->SetLineStyle(2);
   hPt_cohFD__16->SetLineWidth(3);
   hPt_cohFD__16->GetXaxis()->SetTitle("#it{m}_{supcl pair} [GeV/#it{c}^{2}]");
   hPt_cohFD__16->GetXaxis()->SetDecimals();
   hPt_cohFD__16->GetXaxis()->SetLabelFont(42);
   hPt_cohFD__16->GetXaxis()->SetTitleOffset(1.2);
   hPt_cohFD__16->GetXaxis()->SetTitleFont(42);
   hPt_cohFD__16->GetYaxis()->SetTitle("Counts per 80 MeV [-]");
   hPt_cohFD__16->GetYaxis()->SetDecimals();
   hPt_cohFD__16->GetYaxis()->SetNdivisions(3000510);
   hPt_cohFD__16->GetYaxis()->SetLabelFont(42);
   hPt_cohFD__16->GetYaxis()->SetTitleFont(42);
   hPt_cohFD__16->GetZaxis()->SetLabelFont(42);
   hPt_cohFD__16->GetZaxis()->SetTitleOffset(1);
   hPt_cohFD__16->GetZaxis()->SetTitleFont(42);
   hPt_cohFD__16->Draw("HIST SAME");
   
   TH1F *hPt_incFD__17 = new TH1F("hPt_incFD__17","inc #psi' #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow e^{+}e^{-} + #pi^{+}#pi^{-}",40,0,2);
   hPt_incFD__17->SetBinContent(1,116.8019);
   hPt_incFD__17->SetBinContent(2,194.6698);
   hPt_incFD__17->SetBinContent(3,438.007);
   hPt_incFD__17->SetBinContent(4,530.4751);
   hPt_incFD__17->SetBinContent(5,593.7427);
   hPt_incFD__17->SetBinContent(6,764.0788);
   hPt_incFD__17->SetBinContent(7,686.2109);
   hPt_incFD__17->SetBinContent(8,637.5435);
   hPt_incFD__17->SetBinContent(9,754.3453);
   hPt_incFD__17->SetBinContent(10,657.0104);
   hPt_incFD__17->SetBinContent(11,671.6107);
   hPt_incFD__17->SetBinContent(12,452.6072);
   hPt_incFD__17->SetBinContent(13,423.4067);
   hPt_incFD__17->SetBinContent(14,360.139);
   hPt_incFD__17->SetBinContent(15,292.0046);
   hPt_incFD__17->SetBinContent(16,194.6698);
   hPt_incFD__17->SetBinContent(17,165.4693);
   hPt_incFD__17->SetBinContent(18,136.2688);
   hPt_incFD__17->SetBinContent(19,111.9351);
   hPt_incFD__17->SetBinContent(20,53.53418);
   hPt_incFD__17->SetBinContent(21,82.73465);
   hPt_incFD__17->SetBinContent(22,24.33372);
   hPt_incFD__17->SetBinContent(23,14.60023);
   hPt_incFD__17->SetBinContent(24,19.46698);
   hPt_incFD__17->SetBinContent(25,14.60023);
   hPt_incFD__17->SetBinContent(26,4.866744);
   hPt_incFD__17->SetBinContent(28,4.866744);
   hPt_incFD__17->SetBinError(1,23.84208);
   hPt_incFD__17->SetBinError(2,30.77999);
   hPt_incFD__17->SetBinError(3,46.16999);
   hPt_incFD__17->SetBinError(4,50.8103);
   hPt_incFD__17->SetBinError(5,53.75494);
   hPt_incFD__17->SetBinError(6,60.98013);
   hPt_incFD__17->SetBinError(7,57.78938);
   hPt_incFD__17->SetBinError(8,55.70243);
   hPt_incFD__17->SetBinError(9,60.59047);
   hPt_incFD__17->SetBinError(10,56.54645);
   hPt_incFD__17->SetBinError(11,57.1713);
   hPt_incFD__17->SetBinError(12,46.93318);
   hPt_incFD__17->SetBinError(13,45.39397);
   hPt_incFD__17->SetBinError(14,41.86531);
   hPt_incFD__17->SetBinError(15,37.69764);
   hPt_incFD__17->SetBinError(16,30.77999);
   hPt_incFD__17->SetBinError(17,28.37775);
   hPt_incFD__17->SetBinError(18,25.75239);
   hPt_incFD__17->SetBinError(19,23.34008);
   hPt_incFD__17->SetBinError(20,16.14116);
   hPt_incFD__17->SetBinError(21,20.0661);
   hPt_incFD__17->SetBinError(22,10.88237);
   hPt_incFD__17->SetBinError(23,8.429448);
   hPt_incFD__17->SetBinError(24,9.733488);
   hPt_incFD__17->SetBinError(25,8.429448);
   hPt_incFD__17->SetBinError(26,4.866744);
   hPt_incFD__17->SetBinError(28,4.866744);
   hPt_incFD__17->SetMinimum(1);
   hPt_incFD__17->SetMaximum(842.3969);
   hPt_incFD__17->SetEntries(1726);
   hPt_incFD__17->SetStats(0);

   ci = TColor::GetColor("#ff0000");
   hPt_incFD__17->SetLineColor(ci);
   hPt_incFD__17->SetLineStyle(2);
   hPt_incFD__17->SetLineWidth(3);
   hPt_incFD__17->GetXaxis()->SetTitle("#it{m}_{supcl pair} [GeV/#it{c}^{2}]");
   hPt_incFD__17->GetXaxis()->SetDecimals();
   hPt_incFD__17->GetXaxis()->SetLabelFont(42);
   hPt_incFD__17->GetXaxis()->SetTitleOffset(1.2);
   hPt_incFD__17->GetXaxis()->SetTitleFont(42);
   hPt_incFD__17->GetYaxis()->SetTitle("Counts per 80 MeV [-]");
   hPt_incFD__17->GetYaxis()->SetDecimals();
   hPt_incFD__17->GetYaxis()->SetNdivisions(3000510);
   hPt_incFD__17->GetYaxis()->SetLabelFont(42);
   hPt_incFD__17->GetYaxis()->SetTitleFont(42);
   hPt_incFD__17->GetZaxis()->SetLabelFont(42);
   hPt_incFD__17->GetZaxis()->SetTitleOffset(1);
   hPt_incFD__17->GetZaxis()->SetTitleFont(42);
   hPt_incFD__17->Draw("HIST SAME");
   
   TLegend *leg = new TLegend(0.18,0.84,0.6,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.038);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","ALICE Simulation, Pb#minusPb UPC #sqrt{#it{s}_{NN}} = 5.5 TeV, #it{L}_{int} = 7 nb^{-1}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","STARlight, J/#psi and #psi(2#it{S}) #rightarrow e^{+}e^{-}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   leg = new TLegend(0.585,0.575,0.95,0.845,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.036);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("NULL","3.4 < #it{y} < 5.8","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","2.8 < #it{m}_{cl pair} < 3.4 GeV/#it{c}^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPt_all","data","EP");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPt_cohJpsi","coh J/#psi","L");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPt_incJpsi","inc J/#psi","L");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPt_cohFD","feed-down from coh #psi'","L");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPt_incFD","feed-down from inc #psi'","L");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   cPN->Modified();
   cPN->cd();
   cPN->SetSelected(cPN);
}
