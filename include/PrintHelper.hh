// -*- C++ -*-

#ifndef PRINT_HELPER_HH
#define PRINT_HELPER_HH

#include <iomanip>
#include <iostream>
#include <string>

//_____________________________________________________________________________
class PrintHelper
{
public:
  PrintHelper(std::size_t precision,
              std::ios::fmtflags flags,
              std::ostream& ost=std::cout)
    : m_ostream(ost),
      m_precision(ost.precision()),
      m_flags(ost.flags())
  {
    m_ostream.precision(precision);
    m_ostream.setf(flags);
  }
  ~PrintHelper(void)
  {
    m_ostream.precision(m_precision);
    m_ostream.setf(m_flags);
  }

private:
  std::ostream&      m_ostream;
  std::size_t        m_precision;
  std::ios::fmtflags m_flags;

public:
  void Precision(std::size_t precision);
  void Setf(std::ios::fmtflags flags);
  void Set(std::size_t precision, std::ios::fmtflags flags);
};

//_____________________________________________________________________________
inline void
PrintHelper::Precision(std::size_t precision)
{
  m_ostream.precision(precision);
}

//_____________________________________________________________________________
inline void
PrintHelper::Setf(std::ios::fmtflags flags)
{
  m_ostream.setf(flags);
}

//_____________________________________________________________________________
inline void
PrintHelper::Set(std::size_t precision, std::ios::fmtflags flags)
{
  m_ostream.precision(precision);
  m_ostream.setf(flags);
}

#endif
