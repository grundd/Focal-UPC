#ifdef __CLING__
#pragma cling optimize(0)
#endif
void massFit_log()
{
//=========Macro generated from canvas: cPNLog/
//=========  (Sat Aug 26 16:27:19 2023) by ROOT version 6.28/04
   TCanvas *cPNLog = new TCanvas("cPNLog", "",78,53,900,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cPNLog->SetHighLightColor(2);
   cPNLog->Range(0.9333334,0.07057359,4.266666,5.30721);
   cPNLog->SetFillColor(0);
   cPNLog->SetBorderMode(0);
   cPNLog->SetBorderSize(2);
   cPNLog->SetLogy();
   cPNLog->SetLeftMargin(0.14);
   cPNLog->SetRightMargin(0.02);
   cPNLog->SetTopMargin(0.06);
   cPNLog->SetBottomMargin(0.12);
   cPNLog->SetFrameBorderMode(0);
   cPNLog->SetFrameBorderMode(0);
   
   TH1D *frame_600000081bc0__11 = new TH1D("frame_600000081bc0__11","inv mass fit",100,1.4,4.2);
   frame_600000081bc0__11->SetBinContent(1,25069.99);
   frame_600000081bc0__11->SetMinimum(5);
   frame_600000081bc0__11->SetMaximum(98403.88);
   frame_600000081bc0__11->SetEntries(1);
   frame_600000081bc0__11->SetDirectory(nullptr);
   frame_600000081bc0__11->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   frame_600000081bc0__11->SetLineColor(ci);
   frame_600000081bc0__11->GetXaxis()->SetTitle("#it{m}_{cl pair} (GeV/#it{c}^{2})");
   frame_600000081bc0__11->GetXaxis()->SetLabelFont(42);
   frame_600000081bc0__11->GetXaxis()->SetLabelSize(0.05);
   frame_600000081bc0__11->GetXaxis()->SetTitleSize(0.05);
   frame_600000081bc0__11->GetXaxis()->SetTitleOffset(1.05);
   frame_600000081bc0__11->GetXaxis()->SetTitleFont(42);
   frame_600000081bc0__11->GetYaxis()->SetTitle("counts per 80 MeV/#it{c}^{2}");
   frame_600000081bc0__11->GetYaxis()->SetNdivisions(3000510);
   frame_600000081bc0__11->GetYaxis()->SetLabelFont(42);
   frame_600000081bc0__11->GetYaxis()->SetLabelSize(0.05);
   frame_600000081bc0__11->GetYaxis()->SetTitleSize(0.05);
   frame_600000081bc0__11->GetYaxis()->SetTitleOffset(1.25);
   frame_600000081bc0__11->GetYaxis()->SetTitleFont(42);
   frame_600000081bc0__11->GetZaxis()->SetLabelFont(42);
   frame_600000081bc0__11->GetZaxis()->SetTitleOffset(1);
   frame_600000081bc0__11->GetZaxis()->SetTitleFont(42);
   frame_600000081bc0__11->Draw("FUNC");
   
   Double_t fDataHist_fx3006[34] = { 1.441176, 1.523529, 1.605882, 1.688235, 1.770588, 1.852941, 1.935294, 2.017647, 2.1, 2.182353, 2.264706, 2.347059, 2.429412, 2.511765, 2.594118, 2.67647, 2.758823,
   2.841176, 2.923529, 3.005882, 3.088235, 3.170588, 3.252941, 3.335294, 3.417647, 3.5, 3.582353, 3.664706, 3.747059, 3.829412, 3.911765, 3.994117, 4.07647,
   4.158823 };
   Double_t fDataHist_fy3006[34] = { 221.632, 254.9067, 408.5873, 334.6836, 458.0762, 356.2141, 495.9856, 629.312, 749.1611, 817.4679, 960.0189, 1468.96, 1610.476, 2059.918, 2655.933, 4205.045, 6816.271,
   9788.854, 15659.34, 21352.42, 23429.5, 17649.69, 8229.021, 2337.374, 280.8701, 152.8179, 172.1442, 203.3557, 191.3272, 144.2131, 63.1685, 10.20689, 3.052503,
   0 };
   Double_t fDataHist_felx3006[34] = { 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647,
   0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647,
   0.04117647 };
   Double_t fDataHist_fely3006[34] = { 42.66362, 47.09793, 59.43799, 54.06638, 62.71703, 55.50256, 64.57544, 73.72659, 80.06141, 84.33589, 90.76663, 112.9403, 117.6455, 132.334, 150.5628, 189.4108, 241.347,
   288.9885, 365.5356, 426.1566, 446.6892, 387.7886, 265.205, 139.7396, 42.61916, 18.16976, 10.14115, 11.05703, 10.70931, 9.301151, 6.144427, 2.478158, 1.365121,
   0 };
   Double_t fDataHist_fehx3006[34] = { 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647,
   0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647, 0.04117647,
   0.04117647 };
   Double_t fDataHist_fehy3006[34] = { 42.66362, 47.09793, 59.43799, 54.06638, 62.71703, 55.50256, 64.57544, 73.72659, 80.06141, 84.33589, 90.76663, 112.9403, 117.6455, 132.334, 150.5628, 189.4108, 241.347,
   288.9885, 365.5356, 426.1566, 446.6892, 387.7886, 265.205, 139.7396, 42.61916, 18.16976, 10.14115, 11.05703, 10.70931, 9.301151, 6.144427, 2.478158, 1.365121,
   0 };
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(34,fDataHist_fx3006,fDataHist_fy3006,fDataHist_felx3006,fDataHist_fehx3006,fDataHist_fely3006,fDataHist_fehy3006);
   grae->SetName("fDataHist");
   grae->SetTitle("Histogram of fDataHist_plot__fM");
   grae->SetFillStyle(1000);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(24);
   grae->SetMarkerSize(1.2);
   
   TH1F *Graph_fDataHist3006 = new TH1F("Graph_fDataHist3006","Histogram of fDataHist_plot__fM",100,1.12,4.48);
   Graph_fDataHist3006->SetMinimum(0);
   Graph_fDataHist3006->SetMaximum(26263.8);
   Graph_fDataHist3006->SetDirectory(nullptr);
   Graph_fDataHist3006->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_fDataHist3006->SetLineColor(ci);
   Graph_fDataHist3006->GetXaxis()->SetLabelFont(42);
   Graph_fDataHist3006->GetXaxis()->SetTitleOffset(1);
   Graph_fDataHist3006->GetXaxis()->SetTitleFont(42);
   Graph_fDataHist3006->GetYaxis()->SetLabelFont(42);
   Graph_fDataHist3006->GetYaxis()->SetTitleFont(42);
   Graph_fDataHist3006->GetZaxis()->SetLabelFont(42);
   Graph_fDataHist3006->GetZaxis()->SetTitleOffset(1);
   Graph_fDataHist3006->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_fDataHist3006);
   
   grae->Draw("p");
   
   Double_t DSCB_fx8[131] = { 1.371972, 1.372, 1.4, 1.428, 1.456, 1.484, 1.512, 1.54, 1.568, 1.596, 1.624, 1.652, 1.68, 1.708, 1.736, 1.764, 1.792,
   1.82, 1.848, 1.876, 1.904, 1.932, 1.96, 1.988, 2.016, 2.044, 2.072, 2.1, 2.128, 2.156, 2.184, 2.212, 2.24,
   2.268, 2.296, 2.324, 2.352, 2.38, 2.408, 2.436, 2.464, 2.492, 2.52, 2.548, 2.576, 2.604, 2.632, 2.66, 2.688,
   2.716, 2.744, 2.772, 2.786, 2.8, 2.814, 2.828, 2.842, 2.856, 2.87, 2.884, 2.912, 2.94, 2.954, 2.968, 2.982,
   2.996, 3.01, 3.024, 3.038, 3.052, 3.066, 3.08, 3.087, 3.094, 3.101, 3.108, 3.115, 3.122, 3.129, 3.136, 3.143,
   3.15, 3.164, 3.178, 3.192, 3.22, 3.234, 3.248, 3.262, 3.276, 3.29, 3.304, 3.318, 3.332, 3.346, 3.36, 3.374,
   3.388, 3.402, 3.416, 3.444, 3.472, 3.5, 3.528, 3.556, 3.584, 3.612, 3.64, 3.668, 3.696, 3.724, 3.752, 3.78,
   3.808, 3.836, 3.864, 3.892, 3.92, 3.948, 3.976, 4.004, 4.032, 4.06, 4.088, 4.116, 4.144, 4.172, 4.2, 4.2,
   4.228, 4.228028 };
   Double_t DSCB_fy8[131] = { 0, 206.2695, 206.2695, 214.8363, 223.8943, 233.4796, 243.6317, 254.3935, 265.8123, 277.9394, 290.8314, 304.5503, 319.1645, 334.7489, 351.3867, 369.1693, 388.1984,
   408.5864, 430.4585, 453.9539, 479.2281, 506.4548, 535.8287, 567.5686, 601.9207, 639.1629, 679.6096, 723.6177, 771.5931, 823.9992, 881.3663, 944.3036, 1013.513,
   1089.806, 1174.126, 1267.569, 1371.423, 1487.197, 1616.677, 1761.979, 1925.63, 2110.662, 2320.732, 2560.284, 2834.751, 3150.825, 3516.813, 3943.116, 4442.872,
   5032.848, 5734.663, 6576.535, 7061.205, 7595.764, 8186.888, 8842.364, 9571.313, 10384.47, 11294.51, 12316.52, 14658.23, 17020.01, 18152.63, 19228.54, 20229.27,
   21136.89, 21934.56, 22607.05, 23141.2, 23526.37, 23754.78, 23821.09, 23763.86, 23622.68, 23399.08, 23095.39, 22714.8, 22261.26, 21739.4, 21154.49, 20512.31,
   19819.1, 18306.01, 16669.85, 14965.75, 11558.91, 9944.178, 8434.304, 7052.743, 5814.271, 4725.643, 3786.649, 2991.421, 2329.853, 1788.99, 1354.303, 1010.769,
   743.7326, 539.523, 385.8616, 189.1298, 87.57946, 38.31407, 15.83536, 6.183173, 2.280913, 0.794913, 0.2617244, 0.08141093, 0.02392406, 0.006642029, 0.00174213, 0.0004316919,
   0.0001010604, 2.235125e-05, 4.670211e-06, 9.219023e-07, 1.719282e-07, 3.029165e-08, 5.042112e-09, 7.928953e-10, 1.177967e-10, 1.682235e-11, 8.315226e-12, 6.27556e-12, 5.247121e-12, 4.600432e-12, 4.145874e-12, 4.145874e-12,
   4.145874e-12, 0 };
   TGraph *graph = new TGraph(131,DSCB_fx8,DSCB_fy8);
   graph->SetName("DSCB");
   graph->SetTitle("Projection of J/psi DSCB and psi' DSCB");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineStyle(2);
   graph->SetLineWidth(3);
   
   TH1F *Graph_DSCB8 = new TH1F("Graph_DSCB8","Projection of J/psi DSCB and psi' DSCB",131,1.086366,4.513633);
   Graph_DSCB8->SetMinimum(0);
   Graph_DSCB8->SetMaximum(26203.2);
   Graph_DSCB8->SetDirectory(nullptr);
   Graph_DSCB8->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_DSCB8->SetLineColor(ci);
   Graph_DSCB8->GetXaxis()->SetLabelFont(42);
   Graph_DSCB8->GetXaxis()->SetTitleOffset(1);
   Graph_DSCB8->GetXaxis()->SetTitleFont(42);
   Graph_DSCB8->GetYaxis()->SetLabelFont(42);
   Graph_DSCB8->GetYaxis()->SetTitleFont(42);
   Graph_DSCB8->GetZaxis()->SetLabelFont(42);
   Graph_DSCB8->GetZaxis()->SetTitleOffset(1);
   Graph_DSCB8->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_DSCB8);
   
   graph->Draw("l");
   
   Double_t DSCB2_fx9[131] = { 1.371972, 1.372, 1.4, 1.428, 1.456, 1.484, 1.512, 1.54, 1.568, 1.596, 1.624, 1.652, 1.68, 1.708, 1.736, 1.764, 1.792,
   1.82, 1.848, 1.876, 1.904, 1.932, 1.96, 1.988, 2.016, 2.044, 2.072, 2.1, 2.128, 2.156, 2.184, 2.212, 2.24,
   2.268, 2.296, 2.324, 2.352, 2.38, 2.408, 2.436, 2.464, 2.492, 2.52, 2.548, 2.576, 2.604, 2.632, 2.66, 2.688,
   2.716, 2.744, 2.772, 2.8, 2.828, 2.856, 2.884, 2.912, 2.94, 2.968, 2.996, 3.024, 3.052, 3.08, 3.108, 3.136,
   3.164, 3.192, 3.22, 3.248, 3.276, 3.304, 3.332, 3.36, 3.388, 3.416, 3.43, 3.444, 3.458, 3.472, 3.486, 3.5,
   3.514, 3.528, 3.542, 3.556, 3.584, 3.598, 3.612, 3.626, 3.64, 3.654, 3.668, 3.682, 3.696, 3.71, 3.717, 3.724,
   3.731, 3.738, 3.745, 3.752, 3.759, 3.766, 3.773, 3.78, 3.794, 3.808, 3.836, 3.864, 3.878, 3.892, 3.906, 3.92,
   3.934, 3.948, 3.962, 3.976, 3.99, 4.004, 4.018, 4.032, 4.046, 4.06, 4.088, 4.116, 4.144, 4.172, 4.2, 4.2,
   4.228, 4.228028 };
   Double_t DSCB2_fy9[131] = { 0, 0.7318539, 0.7318539, 0.7546955, 0.7785175, 0.8033738, 0.8293218, 0.8564229, 0.8847426, 0.914351, 0.9453232, 0.9777397, 1.011687, 1.047257, 1.084549, 1.123671, 1.164736,
   1.20787, 1.253205, 1.300886, 1.351068, 1.403919, 1.459622, 1.518373, 1.580388, 1.645897, 1.715156, 1.788438, 1.866045, 1.948303, 2.035571, 2.128242, 2.226746,
   2.331554, 2.443185, 2.562212, 2.689264, 2.825039, 2.970307, 3.125925, 3.292841, 3.472115, 3.664926, 3.872592, 4.096592, 4.338585, 4.600437, 4.884255, 5.192424,
   5.527645, 5.892992, 6.291974, 6.728603, 7.207487, 7.733931, 8.31407, 8.955021, 9.665075, 10.45393, 11.33297, 12.31565, 13.41792, 14.65878, 16.06104, 17.65217,
   19.46553, 21.54187, 23.93134, 26.69612, 29.91396, 33.683, 38.12826, 43.41074, 49.74024, 57.39361, 61.82791, 66.74135, 72.20133, 78.28705, 85.09203, 92.72723,
   101.325, 111.044, 122.0568, 133.635, 156.5163, 167.433, 177.7311, 187.2091, 195.6733, 202.9444, 208.8641, 213.3003, 216.1526, 217.3553, 217.2794, 216.5091,
   215.0234, 212.837, 209.9716, 206.4552, 202.3219, 197.6113, 192.3679, 186.6403, 173.9429, 159.9615, 129.9737, 100.1212, 86.1336, 73.1184, 61.2475, 50.62413,
   41.289, 33.2291, 26.38825, 20.67807, 15.98884, 12.1992, 9.184459, 6.823123, 5.001732, 3.617972, 1.818768, 0.8668018, 0.3916449, 0.1677627, 0.06812843, 0.06812843,
   0.06812843, 0 };
   graph = new TGraph(131,DSCB2_fx9,DSCB2_fy9);
   graph->SetName("DSCB2");
   graph->SetTitle("Projection of J/psi DSCB and psi' DSCB");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cc00ff");
   graph->SetLineColor(ci);
   graph->SetLineStyle(2);
   graph->SetLineWidth(3);
   
   TH1F *Graph_DSCB29 = new TH1F("Graph_DSCB29","Projection of J/psi DSCB and psi' DSCB",131,1.086366,4.513633);
   Graph_DSCB29->SetMinimum(0);
   Graph_DSCB29->SetMaximum(239.0909);
   Graph_DSCB29->SetDirectory(nullptr);
   Graph_DSCB29->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_DSCB29->SetLineColor(ci);
   Graph_DSCB29->GetXaxis()->SetLabelFont(42);
   Graph_DSCB29->GetXaxis()->SetTitleOffset(1);
   Graph_DSCB29->GetXaxis()->SetTitleFont(42);
   Graph_DSCB29->GetYaxis()->SetLabelFont(42);
   Graph_DSCB29->GetYaxis()->SetTitleFont(42);
   Graph_DSCB29->GetZaxis()->SetLabelFont(42);
   Graph_DSCB29->GetZaxis()->SetTitleOffset(1);
   Graph_DSCB29->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_DSCB29);
   
   graph->Draw("l");
   
   Double_t CombinedPDF_fx10[131] = { 1.371972, 1.372, 1.4, 1.428, 1.456, 1.484, 1.512, 1.54, 1.568, 1.596, 1.624, 1.652, 1.68, 1.708, 1.736, 1.764, 1.792,
   1.82, 1.848, 1.876, 1.904, 1.932, 1.96, 1.988, 2.016, 2.044, 2.072, 2.1, 2.128, 2.156, 2.184, 2.212, 2.24,
   2.268, 2.296, 2.324, 2.352, 2.38, 2.408, 2.436, 2.464, 2.492, 2.52, 2.548, 2.576, 2.604, 2.632, 2.66, 2.688,
   2.716, 2.744, 2.772, 2.786, 2.8, 2.814, 2.828, 2.842, 2.856, 2.87, 2.884, 2.912, 2.94, 2.954, 2.968, 2.982,
   2.996, 3.01, 3.024, 3.038, 3.052, 3.066, 3.08, 3.087, 3.094, 3.101, 3.108, 3.115, 3.122, 3.129, 3.136, 3.143,
   3.15, 3.164, 3.178, 3.192, 3.22, 3.234, 3.248, 3.262, 3.276, 3.29, 3.304, 3.318, 3.332, 3.346, 3.36, 3.374,
   3.388, 3.402, 3.416, 3.444, 3.472, 3.5, 3.528, 3.556, 3.584, 3.612, 3.64, 3.668, 3.696, 3.724, 3.752, 3.78,
   3.808, 3.836, 3.864, 3.892, 3.92, 3.948, 3.976, 4.004, 4.032, 4.06, 4.088, 4.116, 4.144, 4.172, 4.2, 4.2,
   4.228, 4.228028 };
   Double_t CombinedPDF_fy10[131] = { 0, 207.0013, 207.0013, 215.591, 224.6728, 234.283, 244.461, 255.25, 266.697, 278.8537, 291.7767, 305.5281, 320.1762, 335.7962, 352.4712, 370.293, 389.3631,
   409.7943, 431.7117, 455.2548, 480.5792, 507.8587, 537.2884, 569.087, 603.5011, 640.8088, 681.3248, 725.4061, 773.4592, 825.9475, 883.4019, 946.4318, 1015.74,
   1092.138, 1176.569, 1270.132, 1374.112, 1490.022, 1619.647, 1765.105, 1928.923, 2114.134, 2324.397, 2564.157, 2838.848, 3155.163, 3521.413, 3948, 4448.065,
   5038.375, 5740.556, 6582.827, 7067.71, 7602.493, 8193.851, 8849.572, 9578.778, 10392.2, 11302.53, 12324.83, 14667.18, 17029.67, 18162.68, 19238.99, 20240.15,
   21148.22, 21946.37, 22619.37, 23154.05, 23539.79, 23768.8, 23835.75, 23778.85, 23638.02, 23414.77, 23111.45, 22731.24, 22278.09, 21756.64, 21172.14, 20530.4,
   19837.63, 18325.47, 16690.32, 14987.29, 11582.84, 9969.44, 8461, 7080.985, 5844.185, 4757.366, 3820.332, 3027.233, 2367.981, 1829.643, 1397.714, 1057.199,
   793.4728, 592.9039, 443.2552, 255.8712, 165.8665, 131.0413, 126.8794, 139.8182, 158.7973, 178.526, 195.935, 208.9455, 216.1765, 216.5158, 206.4569, 186.6408,
   159.9616, 129.9738, 100.1212, 73.1184, 50.62413, 33.2291, 20.67807, 12.1992, 6.823123, 3.617972, 1.818768, 0.8668018, 0.3916449, 0.1677627, 0.06812843, 0.06812843,
   0.06812843, 0 };
   graph = new TGraph(131,CombinedPDF_fx10,CombinedPDF_fy10);
   graph->SetName("CombinedPDF");
   graph->SetTitle("Projection of J/psi DSCB and psi' DSCB");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   
   TH1F *Graph_CombinedPDF10 = new TH1F("Graph_CombinedPDF10","Projection of J/psi DSCB and psi' DSCB",131,1.086366,4.513633);
   Graph_CombinedPDF10->SetMinimum(0);
   Graph_CombinedPDF10->SetMaximum(26219.33);
   Graph_CombinedPDF10->SetDirectory(nullptr);
   Graph_CombinedPDF10->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_CombinedPDF10->SetLineColor(ci);
   Graph_CombinedPDF10->GetXaxis()->SetLabelFont(42);
   Graph_CombinedPDF10->GetXaxis()->SetTitleOffset(1);
   Graph_CombinedPDF10->GetXaxis()->SetTitleFont(42);
   Graph_CombinedPDF10->GetYaxis()->SetLabelFont(42);
   Graph_CombinedPDF10->GetYaxis()->SetTitleFont(42);
   Graph_CombinedPDF10->GetZaxis()->SetLabelFont(42);
   Graph_CombinedPDF10->GetZaxis()->SetTitleOffset(1);
   Graph_CombinedPDF10->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_CombinedPDF10);
   
   graph->Draw("l");
   
   TH1D *frame_600000081bc0__12 = new TH1D("frame_600000081bc0__12","inv mass fit",100,1.4,4.2);
   frame_600000081bc0__12->SetBinContent(1,25069.99);
   frame_600000081bc0__12->SetMinimum(5);
   frame_600000081bc0__12->SetMaximum(98403.88);
   frame_600000081bc0__12->SetEntries(1);
   frame_600000081bc0__12->SetDirectory(nullptr);
   frame_600000081bc0__12->SetStats(0);

   ci = TColor::GetColor("#000099");
   frame_600000081bc0__12->SetLineColor(ci);
   frame_600000081bc0__12->GetXaxis()->SetTitle("#it{m}_{cl pair} (GeV/#it{c}^{2})");
   frame_600000081bc0__12->GetXaxis()->SetLabelFont(42);
   frame_600000081bc0__12->GetXaxis()->SetLabelSize(0.05);
   frame_600000081bc0__12->GetXaxis()->SetTitleSize(0.05);
   frame_600000081bc0__12->GetXaxis()->SetTitleOffset(1.05);
   frame_600000081bc0__12->GetXaxis()->SetTitleFont(42);
   frame_600000081bc0__12->GetYaxis()->SetTitle("counts per 80 MeV/#it{c}^{2}");
   frame_600000081bc0__12->GetYaxis()->SetNdivisions(3000510);
   frame_600000081bc0__12->GetYaxis()->SetLabelFont(42);
   frame_600000081bc0__12->GetYaxis()->SetLabelSize(0.05);
   frame_600000081bc0__12->GetYaxis()->SetTitleSize(0.05);
   frame_600000081bc0__12->GetYaxis()->SetTitleOffset(1.25);
   frame_600000081bc0__12->GetYaxis()->SetTitleFont(42);
   frame_600000081bc0__12->GetZaxis()->SetLabelFont(42);
   frame_600000081bc0__12->GetZaxis()->SetTitleOffset(1);
   frame_600000081bc0__12->GetZaxis()->SetTitleFont(42);
   frame_600000081bc0__12->Draw("AXISSAME");
   
   TLegend *leg = new TLegend(0.18,0.75,0.6,0.93,NULL,"brNDC");
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
   entry=leg->AddEntry("NULL","3.4 < #it{y} < 5.8","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","#it{p}_{T} < 0.2 GeV/#it{c}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   leg = new TLegend(0.18,0.2,0.65,0.38,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.038);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("fDataHist","data","EP");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("CombinedPDF","model","L");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("DSCB","J/#psi Crystal Ball","L");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("DSCB2","#psi(2#it{S}) Crystal Ball","L");

   ci = TColor::GetColor("#cc00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   cPNLog->Modified();
   cPNLog->cd();
   cPNLog->SetSelected(cPNLog);
}
