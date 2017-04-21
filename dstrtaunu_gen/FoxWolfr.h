// -*- C++ -*- ----------------------------------------------------------
//
// Package:     Toolbox
// Module:      FoxWolfr
// 
// Description:
//     Calculate Fox-Wolfram moment
//
// Usage:
//
// Author:      Mikihiko Nakao (KEK)
// $Id: FoxWolfr.h 10002 2007-02-26 06:56:17Z katayama $
//
// Revision history
//
// $Log$
// Revision 1.6  2001/12/24 12:03:33  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.5  2001/02/27 17:25:52  nakao
// protection for including twice
//
// Revision 1.4  2000/03/07 11:21:05  katayama
// compatibility with CC5.0
//
// Revision 1.3  1998/07/16 15:53:31  nakao
// compile with egcs-1.0.3a / gcc-2.8.1
//
// Revision 1.1  1997/09/09 03:27:31  nakao
// first version
//
//
// ----------------------------------------------------------------------

#ifndef __FoxWolfr_h__
#define __FoxWolfr_h__

#include "belle.h"
#include "belleCLHEP/Vector/ThreeVector.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class FoxWolfram {
public:
	FoxWolfram();
	~FoxWolfram() {};

	// H(i) = i-th Fox-Wolfram moment
	double H(int i) {
		return (i < 0 || i > 4) ? 0 : sum[i]; }
	// R(i) = i-th normalized Fox-Wolfram moment
	double R(int i) {
		return (i < 0 || i > 4 || sum[0] == 0) ? 0 : sum[i] / sum[0]; }
	void add(const double mag, const double costh);
protected:
	double sum[5];
};

inline FoxWolfram::FoxWolfram()
{
	for (int i=0; i<5; i++)
		sum[i] = 0;
}

inline void FoxWolfram::add(const double mag, const double costh)
{
	double cost2 = costh * costh;

	sum[0] += mag;
	sum[1] += mag * costh;
	sum[2] += mag * (1.5*cost2 - 0.5);
	sum[3] += mag * costh*(2.5*cost2 - 1.5);
	sum[4] += mag * (4.375*cost2*cost2 - 3.75*cost2 + 0.375);
}

template <class Iterator, class Function>
FoxWolfram foxwolfram(Iterator begin, Iterator end, Function func)
{
	FoxWolfram f;

	Iterator p, q;
	for (p = begin; p != end; p++) {
		const Hep3Vector pvec = func(*p);
		double pmag = pvec.mag();
		for (q = p; q != end; q++) {
			const Hep3Vector qvec = func(*q);
			double mag = pmag * qvec.mag();
			double costh = pvec.dot(qvec) / mag;
			if (p != q) mag *= 2;
			f.add(mag, costh);
		}
	}
	return f;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif // __FoxWolfr_h__
