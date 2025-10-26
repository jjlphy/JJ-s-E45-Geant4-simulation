// binom_cp.cpp
// Usage: ./binom_cp <k> <n> [CL]
// - k: successes
// - n: trials
// - CL: confidence level (default 0.95)
// Computes exact (Clopper–Pearson) interval using ROOT's TEfficiency.

#include <iostream>
#include <iomanip>
#include <cmath>
#include "TEfficiency.h"

int main(int argc, char** argv){
    if(argc < 3 || argc > 4){
        std::cerr << "Usage: " << argv[0] << " <k> <n> [CL]\n"
                  << "  k  = number of successes (>=0)\n"
                  << "  n  = number of trials (>=0)\n"
                  << "  CL = confidence level in (0,1), default 0.95\n";
        return 1;
    }

    long long k = std::stoll(argv[1]);
    long long n = std::stoll(argv[2]);
    double CL = (argc == 4) ? std::stod(argv[3]) : 0.95;

    if(n < 0 || k < 0 || k > n || CL <= 0.0 || CL >= 1.0){
        std::cerr << "[Error] Check inputs: require 0 <= k <= n, 0 < CL < 1.\n";
        return 1;
    }
    if(n == 0){
        std::cout << "n = 0 → proportion undefined.\n";
        return 0;
    }

    //const double alpha = 1.0 - CL;
    const double phat  = static_cast<double>(k) / static_cast<double>(n);

   // Exact (Clopper–Pearson) interval via ROOT
   const double lower = TEfficiency::ClopperPearson(static_cast<int>(n), static_cast<int>(k), CL, /*upper=*/false);
   const double upper = TEfficiency::ClopperPearson(static_cast<int>(n), static_cast<int>(k), CL, /*upper=*/true);

    const double err_minus = phat - lower;
    const double err_plus  = upper - phat;

    // (옵션) Wald 표준오차도 참고용으로 출력
    double se_wald = std::sqrt(std::max(0.0, phat*(1.0 - phat) / static_cast<double>(n)));

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Inputs: k=" << k << ", n=" << n << ", CL=" << CL << "\n";
    std::cout << "Proportion p_hat = " << phat
              << "  (" << phat*100.0 << " %)\n";
    std::cout << "Exact (Clopper–Pearson) " << CL*100.0 << "% CI: ["
              << lower << ", " << upper << "]\n";
    std::cout << "Asymmetric errors:  -" << err_minus << "  +" << err_plus << "\n";
    std::cout << "Pretty format:  p = " << phat
              << " _{-" << err_minus << "}^{+" << err_plus << "}\n";
    std::cout << "(Ref) Wald SE ~ " << se_wald << "  (not valid at k=0 or k=n)\n";

    // 퍼센트 표기도 같이
    std::cout << std::setprecision(4);
    std::cout << "Percent CI: [" << lower*100.0 << "%, " << upper*100.0 << "%]\n";
    std::cout << "Percent pretty: " << phat*100.0 << "% _{-" << err_minus*100.0
              << "%}^{+" << err_plus*100.0 << "%}\n";

    return 0;
}
