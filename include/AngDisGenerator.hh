// -*- C++ -*-

/**
 *  for the angular distribution for the E27 experiment
 */

#ifndef ANG_DIS_GENERATOR_HH
#define ANG_DIS_GENERATOR_HH

#include <G4ThreeVector.hh>

//_____________________________________________________________________________
class AngDisGenerator
{
public:
  AngDisGenerator(G4double cost1=1.0, G4double cost2=-1.0);
  virtual ~AngDisGenerator() {};

protected:
  G4double m_cost1;
  G4double m_cost2;

public:
  virtual G4ThreeVector GenerateDirection() const = 0;
  virtual G4double GetDfuncVal(G4double x) const = 0;
};

//_____________________________________________________________________________
class AGSWave : public AngDisGenerator
{
  // D(x)=1/2
public:
  AGSWave(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGSWave() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double) const { return 0.5; }
};

typedef AGSWave AGUniform;

//_____________________________________________________________________________
class AGPWaveFP : public AngDisGenerator
{
  // D(x)=1/2*(1+x)
public:
  AGPWaveFP(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGPWaveFP() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const { return 0.5*(x+1.); }
};

//_____________________________________________________________________________
class AGPWaveBP : public AngDisGenerator
{
  // D(x)=1/2*(1-x)
public:
  AGPWaveBP(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGPWaveBP() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x);}

private:
  inline G4double Dfunc(G4double x) const { return 0.5*(1.-x); }
};

//_____________________________________________________________________________
class AGDWave1 : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGDWave1(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGDWave1() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const { return 0.5*(x*x+1.); }
};

//_____________________________________________________________________________
class AGSigma1385Zero : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigma1385Zero(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGSigma1385Zero() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const
  {
    return
      (6.1 +
	1.2*x +
	(-5.2)/2.*(3*x*x-1))/10.;
  }
};

//_____________________________________________________________________________
class AGSigma1385Plus : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigma1385Plus(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGSigma1385Plus() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const { return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const
  {
    return
      (30.6394+
	(11.7159)*x+
	(-9.52849)/2.*(3*x*x-1)+
	(-13.9936)/2.*(5*x*x*x-3*x)+
	(-7.13241)/8.*(35*x*x*x*x-30*x*x+3)+
	(-3.61133)/8.*(63*x*x*x*x*x-70*x*x*x+15*x)+
	(-2.21394)/16.*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5))/51.;
  }
};

//_____________________________________________________________________________
class AGLambda1405 : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGLambda1405(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGLambda1405() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const
  {
    return
      (1.64+
	1.02*x+
	1.54/2.*(3*x*x-1)+
	0.96/2.*(5*x*x*x-3*x)+
	0.55/8.*(35*x*x*x*x-30*x*x+3)+
	(-0.56)/8.*(63*x*x*x*x*x-70*x*x*x+15*x))/6.;
  }
};

//_____________________________________________________________________________
class AGLambda : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGLambda(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGLambda() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const { return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const
  {
    return
      (13.9+
	9.8*x+
	20.1/2.*(3*x*x-1)-
	1.3/2.*(5*x*x*x-3*x)+
	12.8/8.*(35*x*x*x*x-30*x*x+3))/53.;
  }
};

//_____________________________________________________________________________
class AGSigmaZ : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigmaZ(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGSigmaZ() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const { return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const
  {
    return
      (9.6+
	5.3*x+
	8.0/2.*(3*x*x-1)+
	16.2/2.*(5*x*x*x-3*x)+
	6.7/8.*(35*x*x*x*x-30*x*x+3)+
	10.4/8.*(63*x*x*x*x*x-70*x*x*x+15*x)+
	8.0/16.*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5))/61.;
  }
};

//_____________________________________________________________________________
class AGSigmaP : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigmaP(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGSigmaP() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const
  {
    return
      (37.5+
	7.5*x+
	19.0/2.*(3*x*x-1)+
	23.1/2.*(5*x*x*x-3*x)+
	40.4/8.*(35*x*x*x*x-30*x*x+3)+
	14.6/8.*(63*x*x*x*x*x-70*x*x*x+15*x)+
	31.3/16.*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5))/190.;
  }
};

//_____________________________________________________________________________
class AGPol : public AngDisGenerator
{
public:
  AGPol(G4double cost1=1.0, G4double cost2=-1.0);
  ~AGPol() {};

  G4ThreeVector GenerateDirection() const;
  G4double GetDfuncVal(G4double x) const{ return Dfunc(x); }

private:
  inline G4double Dfunc(G4double x) const { return (1 - 0.4*x)/2.; }
};

#endif
