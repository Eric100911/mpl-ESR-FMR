#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <cmath>

// FMR_BW_1st_Poly_fit.C
void FMR_BW_1st_Poly_fit() {
    // 创建TGraphErrors对象
    const int n = 18;
    double I_m[n] = {1.752, 1.767, 1.776, 1.792, 1.810, 1.819, 1.835, 1.844, 1.850, 1.862, 1.874, 1.883, 1.891, 1.900, 1.915, 1.934, 1.953, 1.975};
    double I_s[n] = {31.0, 30.5, 30.0, 28.5, 27.0, 25.0, 21.0, 17.0, 13.0, 5.0, 2.0, 8.5, 16.5, 25.0, 30.0, 31.0, 32.0, 32.0};
    double B[n]; // 磁感应强度(mT)
    double err_y[n];
    
    // 转换关系: B = 0.1632*I_m + 0.0152
    for (int i = 0; i < n; i++) {
        B[i] = (0.1632 * I_m[i] + 0.0152) * 1000.0; // 转换为mT
        err_y[i] = 1.0; // 假设y误差为1μA
    }

    TGraphErrors *gr = new TGraphErrors(n, B, I_s, nullptr, err_y);
    gr->SetTitle(";B (mT);I_{s} (#muA)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->CenterTitle();
    gr->GetXaxis()->SetTitleOffset(1.1);
    gr->GetYaxis()->SetTitleOffset(1.2);

    // 定义Breit-Wigner拟合函数: f(B) = a0 + a1*B - A / [ (B - mean)^2 + (Γ/2)^2 ]
    // 其中Γ = FWHM
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x - [2]/((x-[3])*(x-[3]) + [4]*[4]/4)", 280.0, 350.0);
    
    // 设置参数初始值和名称
    fitFunc->SetParName(0, "a_{0}");        // 背景常数项
    fitFunc->SetParName(1, "a_{1}");        // 背景线性项
    fitFunc->SetParName(2, "A");            // 峰幅度
    fitFunc->SetParName(3, "#mu");          // 峰中心位置(B)
    fitFunc->SetParName(4, "#Gamma");       // 峰全宽半高(FWHM in B)

    // 设置初始参数值（根据数据特征估计）
    // 中心位置估计: B = 0.1632 * 1.875 + 0.0152 ≈ 0.321 T = 321 mT
    // 幅度估计: 背景约30μA，谷底约2μA，差值为28μA
    // Breit-Wigner峰在中心处的高度为 A/((Γ/2)^2)
    // 设Γ≈10 mT，则高度 = A/(2.5^2) = A/6.25 ≈ 28 → A ≈ 175
    fitFunc->SetParameters(30.0, 0.0, 175.0, 321.0, 10.0); 
    
    // 设置参数范围（物理约束）
    fitFunc->SetParLimits(2, 0, 10000);    // 幅度 > 0
    fitFunc->SetParLimits(4, 1.0, 30.0);   // 宽度 > 0

    // 执行拟合
    gr->Fit(fitFunc, "R"); // "R" 表示在函数范围内拟合

    // 创建背景函数
    TF1 *bgFunc = new TF1("bgFunc", "[0] + [1]*x", 280.0, 350.0);
    bgFunc->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1));
    bgFunc->SetLineColor(kGreen);
    bgFunc->SetLineStyle(2);
    bgFunc->SetLineWidth(2);

    // 创建画布并绘制结果
    TCanvas *c1 = new TCanvas("c1", "Breit-Wigner Peak with Linear Background", 800, 600);
    
    gr->Draw("AP");
    fitFunc->Draw("same");
    bgFunc->Draw("same");
    
    // 添加图例
    TLegend *leg = new TLegend(0.7, 0.50, 0.9, 0.70);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->AddEntry(gr, "Data", "lep");
    leg->AddEntry(fitFunc, "Total Fit", "l");
    leg->AddEntry(bgFunc, "Background", "l");
    leg->Draw();

    // 显示拟合参数
    TPaveText *pt = new TPaveText(0.13, 0.13, 0.42, 0.37, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->AddText(Form("#it{#mu} = %.2f #pm %.2f mT", fitFunc->GetParameter(3), fitFunc->GetParError(3)));
    pt->AddText(Form("#it{#Gamma} = %.2f #pm %.2f mT", fitFunc->GetParameter(4), fitFunc->GetParError(4)));
    pt->AddText(Form("#chi^{2}/NDF = %.2f", fitFunc->GetChisquare()/fitFunc->GetNDF()));
    pt->Draw();

    // 保存图片
    c1->SaveAs("FMR_BreitWigner_1st_Fit.png");
}