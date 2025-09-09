// check_segments.C
// 사용법:
// root -l
// root [0] .L gen.dict.C+
// root [1] .L check_segments.C+
// root [2] check_segments("E45_BVH2.root");

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include <vector>
#include <iostream>
#include <algorithm>

#define SEG(P) (P.GetMother(0))  // 세그 ID는 Mother(0)로 확정

void check_segments(const char* filename="E45_BVH2.root",
                    int nU=26, int nD=40, bool verbose=true)
{
    TFile* f = TFile::Open(filename,"READ");
    if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
    TTree* tr = (TTree*)f->Get("g4hyptpc");
    if(!tr){ std::cerr<<"[ERR] no TTree g4hyptpc\n"; return; }

    std::vector<TParticle> *bvhu=nullptr, *bvhd=nullptr;
    tr->SetBranchAddress("BVH_U",&bvhu);
    tr->SetBranchAddress("BVH_D",&bvhd);

    std::vector<long long> cntU(nU,0), cntD(nD,0);
    int minU=+1e9, maxU=-1e9, minD=+1e9, maxD=-1e9;

    const Long64_t N = tr->GetEntries();
    for(Long64_t i=0;i<N;++i){
        tr->GetEntry(i);
        if(bvhu) for(const auto& p:*bvhu){
            int id=SEG(p);
            minU=std::min(minU,id); maxU=std::max(maxU,id);
            if(0<=id && id<nU) cntU[id]++;
        }
        if(bvhd) for(const auto& p:*bvhd){
            int id=SEG(p);
            minD=std::min(minD,id); maxD=std::max(maxD,id);
            if(0<=id && id<nD) cntD[id]++;
        }
    }

    if(minU>maxU){minU=0; maxU=-1;}
    if(minD>maxD){minD=0; maxD=-1;}

    auto print_cov = [&](const char* name, const std::vector<long long>& cnt, int n, int minId, int maxId){
        int seen=0; for(int i=0;i<n;++i) if(cnt[i]>0) ++seen;
        double cov = 100.0*seen/std::max(1,n);
        std::cout << "\n["<<name<<"]\n";
        std::cout << "  Defined range: 0.."<<(n-1)<<"\n";
        if(minId<=maxId) std::cout << "  Observed ID range: "<<minId<<".."<<maxId<<"\n";
        else              std::cout << "  Observed ID range: (none)\n";
        std::cout << "  Seen segments: "<<seen<<"/"<<n<<" ("<<cov<<"%)\n";
        std::cout << "  Missing IDs  : ";
        bool first=true;
        for(int i=0;i<n;++i) if(cnt[i]==0){ if(!first) std::cout<<","; std::cout<<i; first=false; }
        if(first) std::cout << "(none)";
        std::cout << "\n";
        if(verbose){
            std::cout << "  Nonzero counts (id:count): ";
            first=true;
            for(int i=0;i<n;++i) if(cnt[i]>0){
                if(!first) std::cout<<", ";
                std::cout<<i<<":"<<cnt[i];
                first=false;
            }
            if(first) std::cout << "(none)";
            std::cout << "\n";
        }
    };

    print_cov("BVH_U", cntU, nU, minU, maxU);
    print_cov("BVH_D", cntD, nD, minD, maxD);
}
