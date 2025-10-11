// BH2_SCH_CutEfficiency_DEBUG.C
// - Denominator(분모):
//   [BH2 hit] && [HTOF(any) >= 1] && [HTOF event-level exclude: 배제세그 1개라도 포함 시 drop]
//   && [HTOF multiplicity (unique, 배제없이) >= 2]
// - Input cut files: 각 줄에 (BH2, SCH) 한 쌍(다양한 포맷 허용)
//   권장 표기: "(BH2,SCH)=(h, s)"; "h,s" 또는 "h s" 도 허용
// - Diagnostics: CWD, 파일 stat(FileStat_t), 첫 10줄 프리뷰, 범위 밖 좌표 경고(최대 20건)
// - Output: tight/fit/교집합/합집합 비율
//
// 사용 예:
//   root -l
//   .L BH2_SCH_CutEfficiency_DEBUG.C+
//   BH2_SCH_CutEfficiency_DEBUG("E45_piplusn.root","g4hyptpc","15-20", "./BH2&&SCH_tight_cut.txt","./BH2&&SCH_fit_cut.txt",15,64);
//                               "./BH2&&SCH_tight_cut.txt","./BH2&&SCH_fit_cut.txt",
//                               15,64);

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TSystemFile.h>

#include <unordered_set>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <cmath>

// ===== HTOF geometry =====
static const int    kNPlaneHTOF  = 8;
static const int    kSegPerPlane = 4;
static const int    kNSegHTOF    = kNPlaneHTOF * kSegPerPlane; // 32
static const double kPitch       = 68.0;   // [mm]
static const double kL           = 337.0;  // [mm]
static const double kHTOFx       = 0.0;
static const double kHTOFy       = 12.0;
static const double kHTOFz       = 0.0;

struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

// -------- HTOF mapping helpers --------
inline void toLocalCenter(double x,double y,double z,double& lx,double& ly,double& lz){
  lx = x - kHTOFx;  ly = y - kHTOFy;  lz = z - kHTOFz;
}
inline void planeAxes(int i, double& exx, double& exz, double& ezx, double& ezz){
  const double th = i * 2.0 * TMath::Pi() / kNPlaneHTOF; // i*45°
  exx =  std::cos(th);  exz = -std::sin(th);
  ezx =  std::sin(th);  ezz =  std::cos(th);
}
inline int inferLocalFromX(double xloc, bool flip_local=false){
  double idx_f = xloc / kPitch + 1.5;
  int j = TMath::Nint(idx_f);
  if(j<0) j=0; if(j>=kSegPerPlane) j=kSegPerPlane-1;
  if(flip_local) j = (kSegPerPlane-1)-j;
  return j;
}
inline bool mapHitToHTOFseg(double x,double y,double z,
                            int L_sign,int plane_offset,bool flip_local,
                            int& seg,double R_tol_mm)
{
  double lx,ly,lz; toLocalCenter(x,y,z,lx,ly,lz);
  const double R = std::hypot(lx,lz);
  if(R_tol_mm>0 && std::abs(R - kL) > R_tol_mm) return false;

  int best_i=-1; double best_abs=1e99, best_xloc=0.0;
  for(int i=0;i<kNPlaneHTOF;++i){
    double exx,exz,ezx,ezz; planeAxes(i,exx,exz,ezx,ezz);
    const double zloc = lx*ezx + lz*ezz;
    const double diff = std::abs(zloc - (L_sign*kL));
    if(diff < best_abs){ best_abs = diff; best_i = i; best_xloc = lx*exx + lz*exz; }
  }
  if(best_i<0) return false;

  int plane = (best_i + plane_offset) % kNPlaneHTOF;
  if(plane<0) plane += kNPlaneHTOF;
  int local = inferLocalFromX(best_xloc, flip_local);
  seg   = plane*kSegPerPlane + local;
  return (seg>=0 && seg<kNSegHTOF);
}

// -------- spec parser: "15-20,22" / "none" --------
static std::unordered_set<int> parse_spec(const char* spec, int lo, int hi){
  std::unordered_set<int> out;
  if(!spec) return out;
  std::string s(spec);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  if(s.empty() || s=="none" || s=="NONE") return out;
  std::stringstream ss(s); std::string tok;
  while(std::getline(ss, tok, ',')){
    if(tok.empty()) continue;
    auto pos = tok.find('-');
    if(pos==std::string::npos){
      int v = std::stoi(tok); if(v>=lo && v<=hi) out.insert(v);
    }else{
      int a = std::stoi(tok.substr(0,pos));
      int b = std::stoi(tok.substr(pos+1));
      if(a>b) std::swap(a,b);
      a = std::max(a, lo); b = std::min(b, hi);
      for(int v=a; v<=b; ++v) out.insert(v);
    }
  }
  return out;
}

// -------- unique segments from TParticle branch --------
static inline void unique_seg_from_branch(const std::vector<TParticle>* v,
                                          int nmax, double ecutMeV,
                                          std::set<int>& out)
{
  out.clear();
  if(!v) return;
  for(const auto& p : *v){
    if(ecutMeV>0.0 && p.GetWeight()<ecutMeV) continue;
    int seg = p.GetMother(1);
    if(0<=seg && seg<nmax) out.insert(seg);
  }
}

// -------- file diagnostics using FileStat_t (portable) --------
static void print_file_head(const char* path, int nlines=10){
  FileStat_t st;
  int rc = gSystem->GetPathInfo(path, st);
  if(rc==0) std::cout << "[FILE] " << path << "  size="<< st.fSize <<" bytes\n";
  else      std::cout << "[FILE] " << path << "  (stat failed)\n";

  std::ifstream fin(path);
  if(!fin.is_open()){
    std::cout << "[FILE] cannot open\n";
    return;
  }
  std::cout << "[HEAD] first " << nlines << " line(s):\n";
  std::string line; int c=0;
  while(c<nlines && std::getline(fin,line)){
    std::cout << "  " << line << "\n"; c++;
  }
  fin.close();
}

// -------- robust cut loader (FIXED) --------
// - 등호(=) 오른쪽의 마지막 두 정수만 사용하여 레이블 숫자 오염 방지
// - "(BH2,SCH)=(h, s)", "h,s", "h s" 등 다양한 포맷 허용
static bool load_cut_pairs(const char* path, int BH2_NSEG, int SCH_NSEG,
                           std::unordered_set<long long>& pairs, int& n_lines_ok)
{
  pairs.clear(); n_lines_ok=0;
  std::ifstream fin(path);
  if(!fin.is_open()){
    std::cerr << "[ERR] cannot open cut file: " << path << "\n";
    return false;
  }
  std::string line;
  static int warnOutRange=0; const int warnCap=20;

  while(std::getline(fin, line)){
    if(line.empty()) continue;

    // CR 제거 + 앞뒤 공백 제거
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    auto p0 = line.find_first_not_of(" \t");
    if(p0==std::string::npos) continue;
    if(line[p0]=='#' || line[p0]=='%') continue;

    int h=-1, s=-1;
    bool ok=false;

    // 1) 정식 포맷 "(BH2,SCH)=(h, s)" 시도 (공백/콤마 허용)
    if(sscanf(line.c_str(), " (BH2 , SCH ) = ( %d , %d )", &h, &s)==2 ||
       sscanf(line.c_str(), " (BH2,SCH)=(%d,%d)", &h, &s)==2 ||
       sscanf(line.c_str(), " (BH2,SCH) = ( %d , %d )", &h, &s)==2) {
      ok = true;
    }

    // 2) 콤마/공백 두 숫자만 있는 라인 허용: "h,s" 또는 "h s"
    if(!ok){
      // 등호 뒤만 사용 (레이블의 숫자 무시)
      const char* rhs_c = strchr(line.c_str(), '=');
      std::string rhs = rhs_c ? std::string(rhs_c+1) : line;

      // 숫자/공백/콤마만 남기기
      std::string buf; buf.reserve(rhs.size());
      for(char c: rhs){
        if(std::isdigit((unsigned char)c) || c=='-' || c==' ' || c==',' || c=='\t')
          buf.push_back(c);
      }
      // 마지막 두 정수 취함
      std::vector<long long> nums; nums.reserve(4);
      std::stringstream ss(buf);
      long long v; char ch;
      while(ss >> v){
        nums.push_back(v);
        if(ss.peek()==',') ss >> ch;
      }
      if(nums.size()>=2){
        h = (int)nums[nums.size()-2];
        s = (int)nums[nums.size()-1];
        ok = true;
      }
    }

    if(!ok) continue;

    if(!(0<=h && h<BH2_NSEG && 0<=s && s<SCH_NSEG)){
      if(warnOutRange < warnCap){
        std::cerr << "[WARN] out-of-range pair ignored: (BH2,SCH)=("<<h<<","<<s<<")\n";
        warnOutRange++;
        if(warnOutRange==warnCap) std::cerr << "[WARN] (further out-of-range warnings suppressed)\n";
      }
      continue;
    }

    long long key = ( ( (long long)h & 0xffffffffLL )<<32 ) | ( (long long)s & 0xffffffffLL );
    pairs.insert(key);
    n_lines_ok++;
  }
  fin.close();
  return true;
}

// ================== MAIN (DEBUG) ==================
void BH2_SCH_CutEfficiency_DEBUG(const char* filename="E45_piplusn.root",
                           const char* treename="g4hyptpc",
                           const char* htof_exclude_spec="15-20",
                           const char* cutfile_tight="BH2&&SCH_tight_cut.txt",
                           const char* cutfile_fit  ="BH2&&SCH_fit_cut.txt",
                           int BH2_NSEG=15, int SCH_NSEG=64,
                           double ecutBH2MeV=0.0, double ecutHTOFMeV=0.0, double ecutSCHMeV=0.0,
                           double R_tol_mm=200.0, int L_sign=+1, int plane_offset=0, bool flip_local=false)
{
  // --- CWD & file heads ---
  std::cout << "[CWD] " << gSystem->pwd() << "\n";
  print_file_head(cutfile_tight, 10);
  print_file_head(cutfile_fit,   10);

  // --- load cut pairs ---
  std::unordered_set<long long> CUT_TIGHT, CUT_FIT;
  int n_ok_tight=0, n_ok_fit=0;
  bool okT = load_cut_pairs(cutfile_tight, BH2_NSEG, SCH_NSEG, CUT_TIGHT, n_ok_tight);
  bool okF = load_cut_pairs(cutfile_fit,   BH2_NSEG, SCH_NSEG, CUT_FIT,   n_ok_fit);
  std::cout << "[INFO] loaded from '"<<cutfile_tight<<"': "<<n_ok_tight<<" pairs\n";
  std::cout << "[INFO] loaded from '"<<cutfile_fit  <<"': "<<n_ok_fit  <<" pairs\n";

  // sample print
  auto print_sample = [](const char* tag, const std::unordered_set<long long>& S){
    std::cout << "[SAMPLE] " << tag << " (up to 20 pairs):\n";
    int cnt=0;
    for(auto key : S){
      int h = (int)(key>>32);
      int s = (int)(key & 0xffffffffLL);
      std::cout << "  (BH2,SCH)=("<<h<<","<<s<<")\n";
      if(++cnt>=20) break;
    }
    if(S.empty()) std::cout << "  (none)\n";
  };
  print_sample("tight", CUT_TIGHT);
  print_sample("fit  ", CUT_FIT);

  if(!okT || !okF){
    std::cerr << "[ERR] One or both cut files could not be opened. Aborting.\n";
    return;
  }

  // --- ROOT I/O ---
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] Cannot open "<<filename<<"\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] Tree "<<treename<<" not found\n"; f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *HTOF=nullptr, *SCH=nullptr;
  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);
  if(T->GetBranch("SCH")) T->SetBranchAddress("SCH",&SCH);
  else { std::cerr<<"[ERR] Branch 'SCH' not found\n"; f->Close(); return; }

  const auto HTOF_EXCL = parse_spec(htof_exclude_spec, 0, kNSegHTOF-1);

  const Long64_t N_total = T->GetEntries();
  Long64_t N_bh2_any=0, N_pass1=0, N_denom=0;
  Long64_t N_hit_tight=0, N_hit_fit=0, N_hit_both=0;

  // --- event loop ---
  for(Long64_t ie=0; ie<N_total; ++ie){
    T->GetEntry(ie);
    if(!BH2||!HTOF||!SCH) continue;

    // BH2 unique
    std::set<int> bh2_segs;
    unique_seg_from_branch(BH2, BH2_NSEG, ecutBH2MeV, bh2_segs);
    if(bh2_segs.empty()) continue;
    N_bh2_any++;

    // HTOF unique (all, for exclude + multiplicity)
    std::set<int> htof_all;
    for(const auto& p : *HTOF){
      if(ecutHTOFMeV>0.0 && p.GetWeight()<ecutHTOFMeV) continue;
      int seg=-1;
      if(!mapHitToHTOFseg(p.Vx(),p.Vy(),p.Vz(), L_sign,plane_offset,flip_local, seg, R_tol_mm)) continue;
      if(seg<0 || seg>=kNSegHTOF) continue;
      htof_all.insert(seg);
    }
    if(htof_all.empty()) continue;

    // event-level exclude
    bool has_excl=false;
    for(int s: htof_all){ if(HTOF_EXCL.count(s)){ has_excl=true; break; } }
    if(has_excl) continue;
    N_pass1++;

    // denominator: mult >= 2
    if((int)htof_all.size() < 2) continue;
    N_denom++;

    // SCH unique
    std::set<int> sch_unique;
    unique_seg_from_branch(SCH, SCH_NSEG, ecutSCHMeV, sch_unique);
    if(sch_unique.empty()) continue;

    // check cuts
    bool match_tight=false, match_fit=false;
    for(int h : bh2_segs){
      for(int s : sch_unique){
        long long key = ( ((long long)h)<<32 ) | (unsigned long long)s;
        if(!match_tight && CUT_TIGHT.count(key)) match_tight=true;
        if(!match_fit   && CUT_FIT  .count(key)) match_fit=true;
        if(match_tight && match_fit) break;
      }
      if(match_tight && match_fit) break;
    }
    if(match_tight) ++N_hit_tight;
    if(match_fit)   ++N_hit_fit;
    if(match_tight && match_fit) ++N_hit_both;
  }

  auto pct=[](Long64_t a,Long64_t b){ return (b>0)?100.0*double(a)/double(b):0.0; };

  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"\n==================== CUT EFFICIENCY (DEBUG) ====================\n";
  std::cout<<"Input file                         : "<<filename<<"\n";
  std::cout<<"Total events                       : "<<N_total<<"\n";
  std::cout<<"BH2 hit (any seg)                  : "<<N_bh2_any<<"\n";
  std::cout<<"BH2 hit & HTOF>=1 & !{excl}        : "<<N_pass1<<" ("<<pct(N_pass1,N_bh2_any)<<" %)\n";
  std::cout<<"Denominator (HTOF mult>=2, no excl): "<<N_denom<<" ("<<pct(N_denom,N_pass1)<<" %)\n";
  std::cout<<"HTOF excluded segs (event-level)   : {"<<(htof_exclude_spec?htof_exclude_spec:"")<<"}\n";
  std::cout<<"Loaded tight pairs                 : "<<n_ok_tight<<"\n";
  std::cout<<"Loaded fit   pairs                 : "<<n_ok_fit  <<"\n";
  std::cout<<"---------------------------------------------------------------\n";
  std::cout<<"TIGHT  : "<<N_hit_tight<<" / "<<N_denom<<"  = "<<pct(N_hit_tight,N_denom)<<" %\n";
  std::cout<<"FIT    : "<<N_hit_fit  <<" / "<<N_denom<<"  = "<<pct(N_hit_fit  ,N_denom)<<" %\n";
  Long64_t N_both  = N_hit_both;
  Long64_t N_union = N_hit_tight + N_hit_fit - N_both;
  std::cout<<"BOTH(∩): "<<N_both <<" / "<<N_denom<<"  = "<<pct(N_both ,N_denom)<<" %\n";
  std::cout<<"UNION(∪): "<<N_union<<" / "<<N_denom<<"  = "<<pct(N_union,N_denom)<<" %\n";
  std::cout<<"NONE   : "<<(N_denom - N_union)<<" / "<<N_denom<<"  = "<<pct(N_denom - N_union,N_denom)<<" %\n";
  std::cout<<"===============================================================\n";

  f->Close();
}
