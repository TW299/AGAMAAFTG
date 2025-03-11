#pragma once
#include "actions_base.h"
#include "actions_staeckel.h"
#include "actions_newgenfnc.h"
#include "actions_newisochrone.h"
#include "actions_focal_distance_finder.h"
#include "potential_base.h"
#include "potential_utils.h"
#include "math_core.h"
#include "math_specfunc.h"
#include "math_fit.h"
#include "math_linalg.h"
#include "math_spline.h"
#include "math_fourier.h"
#include "coord.h"
#include "particles_base.h"
#include <iostream>
#include <map>
#include <complex>
#include <cassert>

namespace actions {
	/*
	 Maps sphere to prolate ellipsoid extended along z
	*/
	class EXP PointTrans {
	private:
		coord::UVSph cs;
	public:
		PointTrans(double _D) : cs(_D) {}
		PointTrans(coord::UVSph _cs) : cs(_cs) {}
		coord::PosMomSph Cyl2Sph(const coord::PosMomCyl Rz) const;
		coord::PosMomCyl Sph2Cyl(const coord::PosMomSph rp) const;
	};
	/*
	 Combination of a PointTrans and an Isochrone aa map
	*/
	class EXP ToyMap {
	private:
		PointTrans PT;
		Isochrone Is;
	public:
		ToyMap(double _D, double _Js, double _b) : PT(_D), Is(_Js, _b) {}
		ToyMap(coord::UVSph _cs, Isochrone _Is) : PT(_cs), Is(_Is) {}
		Actions pq2J(const coord::PosMomCyl Rzp) const {
			return Is.pq2J(PT.Cyl2Sph(Rzp));
		}
		ActionAngles pq2aa(const coord::PosMomCyl& Rz) const {
			coord::PosMomSph rp(PT.Cyl2Sph(Rz));
			if (Is.H(rp) < 0) return Is.pq2aa(rp);
			printf("ToyMap::pq2aa: H>=0!");
			exit(0);
		}
		coord::PosMomCyl from_aaT(ActionAngles aaT) const {
			return PT.Sph2Cyl(Is.aa2pq(aaT));
		}
	};

	/*
	  Map of equatorial plane of (r,theta,phi) to system in (x,y) plane
	  extended along y
	*/
	class PsiCoord {
	private:
		const double a, b;
	public:
		PsiCoord(const double _a, const double _b) : a(_a), b(_b) {}
		double Psi(const double phi, double& dpsidphi) const {
			dpsidphi = 1 + 2 * a * cos(2 * phi) + 4 * b * cos(4 * phi);
			return phi + a * sin(2 * phi) + b * sin(4 * phi);
		}
		double Phi(const double psi, double& dpsidphi) const {
			double phi = psi, D = Psi(phi, dpsidphi) - psi;
			while (fabs(D) > 1e-5) {
				phi -= D / dpsidphi;
				D = Psi(phi, dpsidphi) - psi;
			}
			return phi;
		}
	};
	class EXP xyPointTrans {
		const double Delta, Delta2;
		PsiCoord PC;
	public:
		xyPointTrans(double D, double a, double b) : Delta(D), Delta2(D* D), PC(a, b) {}
		coord::PosMomCar rp2xp(const coord::PosMomSph rp) const;
		coord::PosMomSph xp2rp(const coord::PosMomCar xp) const;
	};
	class EXP xyToyMap {
	private:
	public:
		Isochrone Is;
		xyPointTrans PT;
		xyToyMap(double _D, double _aPT, double _bPT, double _Js, double _bIso) :
			PT(_D, _aPT, _bPT), Is(_Js, _bIso) {
		}
		double H(const coord::PosMomCar xp) const {
			coord::PosMomSph rp(PT.xp2rp(xp));
			return Is.H(rp);
		}
		Actions pq2J(const coord::PosMomCar xp) const {
			coord::PosMomSph rp(PT.xp2rp(xp));
			if (Is.H(rp) < 0) return Is.pq2J(rp);
			printf("ToyMap::pq2J: H>=0!");
			exit(0);
		}
		ActionAngles pq2aa(const coord::PosMomCar& xp) const {
			coord::PosMomSph rp(PT.xp2rp(xp));
			if (Is.H(rp) < 0) return Is.pq2aa(rp);
			printf("ToyMap::pq2aa: H>=0!");
			exit(0);
		}
		coord::PosMomCar aa2pq(ActionAngles aaT) const {
			return PT.rp2xp(Is.aa2pq(aaT));
		}
		coord::PosMomCar aa2pq(Actions J, Angles thetas) const {
			return aa2pq(ActionAngles(J, thetas));
		}
	};
	/* Class to hold the Fourier decomposition of the residual Hamiltonian
	 * after torus fitting
	*/
	class EXP PerturbingHamiltonian {
	private:
		GenFncIndices indices;
		std::vector<std::complex<double> > values;
	public:
		PerturbingHamiltonian(const GenFncIndices& _indices,
			const std::vector<std::complex<double> >& _values) :
			indices(_indices), values(_values) {
		}
		GenFncIndex index(const int i) const {
			return indices[i];
		}
		std::complex<double> value(const int i) const {
			return values[i];
		}
		int numTerms() {
			return indices.size();
		}
		std::vector<std::complex<double> > get_hn(const GenFncIndex&, std::vector<float>&) const;
		PerturbingHamiltonian& operator *= (const double a);
		PerturbingHamiltonian  operator * (const double a) const;
		PerturbingHamiltonian& operator += (const PerturbingHamiltonian&);
		PerturbingHamiltonian  operator + (const PerturbingHamiltonian&) const;
	};
	/*
	 * Base of all torus classes
	 */
	class EXP Torus {
	public:
		Actions J;
		Frequencies freqs;
		GenFnc GF;
		double E;
		coord::UVSph cs;
		Isochrone Is;
		ToyMap TM;
		GenFncIndices tried;
		/* Creator called by TorusGenerator rather than users */
		Torus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			const Isochrone& _Is, const coord::UVSph& _cs, double _E) :
			J(_J), freqs(_freqs), GF(_GF), Is(_Is), cs(_cs), E(_E), TM(_cs, _Is) {
		}
		/*		Torus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
					  const Isochrone& _Is, const coord::UVSph& _cs, double _E, GenFncIndices _tried) :
					J(_J), freqs(_freqs), GF(_GF), Is(_Is), cs(_cs),
					E(_E), tried(_tried), TM(_cs,_Is) {}*/
		void printGF(void) {
			GF.print();
		}
		coord::PosMomCyl from_toy(const Angles&) const;
		coord::PosMomCyl from_true(const Angles&) const;
		coord::PosMomCyl from_aaT(const ActionAngles&) const;//toy inputs
		coord::PosCyl PosDerivJ(const Angles&, actions::DerivAct<coord::Cyl>&) const;
		coord::PosCyl PosDerivs(const Angles&, actions::DerivAng<coord::Cyl>&,
			double* det = NULL) const;/*dR/dtheta etx*/
		/* compute the surface of section z=0 pz>0 */
		void zSoS(std::vector<double>& R, std::vector<double>& vR, const int N,
			double& Rmin, double& Rmax, double& Vmax, const double z0 = 0) const;
		/* compute the surface of section R=Rbar pR>0 */
		void rSoS(std::vector<double>& z, std::vector<double>& vz, const double Rbar, const int N,
			double& zmax, double& Vmax) const;
		Frequencies Omega(void) const {
			return freqs;
		}
		/* Compute the orbit from given angles using computed frequencies */
		std::vector<std::pair<coord::PosVelCyl, double> > orbit(const Angles& theta0, double dt, double T) const;
		/* Does the torus pass through given point? If so for
		 * what true thetas? */
		bool containsPoint(const coord::PosCyl& pt, std::vector<Angles>& As,
			std::vector<coord::VelCyl>& Vs,
			std::vector<double>& Jacobs,
			std::vector<actions::DerivAngCyl>* dA = NULL,
			const double tol = 1e-6) const;
		/* returns density contributed at location */
		double density(const coord::PosCyl&) const;
		void write(FILE*) const;
		void read(FILE*);
		Torus& operator *= (const double);
		Torus& operator += (const Torus&);
		const Torus operator * (const double) const;
		const Torus operator + (const Torus&) const;
	};

	//interpolate between 2 tori
	EXP Torus interpTorus(const double x, const Torus& T0, const Torus& T1);

	//interpolae on an indexed array of tori 
	EXP Torus interpTorus(const double x, std::vector<double>&, std::vector<Torus>&);

	/* The class of tori that include perturbing Hamiltonians. Can be
	 * interpolated */
	class EXP eTorus : public Torus {
	private:
	public:
		PerturbingHamiltonian pH;
		eTorus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			const Isochrone& _Is, const coord::UVSph& _cs, double _E,
			const PerturbingHamiltonian& _pH) :
			Torus(_J, _freqs, _GF, _Is, _cs, _E), pH(_pH) {
			//printf("eTorus created: %d terms in pH\n", pH.numTerms());
		}
		eTorus(const Torus& T, const PerturbingHamiltonian& _pH) :
			Torus(T), pH(_pH) {
			//printf("eTorus at E = %f created: %d terms in pH\n",T.E, pH.numTerms());
		}
		PerturbingHamiltonian Hns() const {
			return pH;
		}
		std::vector<std::complex<double> > get_hn(const GenFncIndex& Indx, std::vector<float>& multiples) const {
			return pH.get_hn(Indx, multiples);
		}
		eTorus& operator *= (const double);
		eTorus& operator += (const eTorus&);
		const eTorus operator * (const double) const;
		const eTorus operator + (const eTorus&) const;
	};

	//interpolate between 2 tori
	EXP eTorus interpeTorus(const double x, const eTorus& T0, const eTorus& T1);

	//interpolae on an indexed array of tori 
	EXP eTorus interpeTorus(const double x, std::vector<double>&, std::vector<eTorus>&);

	/*
	 * Class for fitting torus to an orbit
	*/
	class EXP TMfitter : public math::IFunctionNoDeriv {
	private:
		double pphi, Delta2, xmin, ymax, xbar, Frat, aPT, bPT;
		std::vector<std::pair<coord::PosMomCyl, double> >& traj;
	public:
		TMfitter(const potential::BasePotential&,
			std::vector<std::pair<coord::PosMomCyl, double> >&, double);
		std::vector<double> fitTM() const;
		virtual double value(double) const;
		virtual unsigned int numVars() const { return 4; }
		virtual unsigned int numValues() const { return 1; }
	};
	/*
	 * Class for generators of tori.
	*/
	class EXP TorusGenerator {
	private:
		const potential::BasePotential& pot;
		const double defaultTol, invPhi0;
		math::LinearInterpolator2d interpD; //for Delta values
		math::LinearInterpolator2d interpR; //for Rshell values
#ifdef TEST
		/* Test_it compares analytic and numerical derivatives */
		void test_it(actions::Actions& J, std::vector<double>&);
#endif
		/* Hamilton computes residual H for Phi+ePhi where
		 * ePhi is an optional additional potential not used in fitting */
		double Hamilton(const Torus&, const potential::BasePotential*, const Angles&);
		/* PerturbingHamiltonian computes Fourier decomposition
		 * of the residual H */
		PerturbingHamiltonian get_pH(const Torus&,
			int nf, bool ifp, const potential::BasePotential*);
		void setConsts(actions::Actions, double&, double&, double&, Isochrone&, coord::UVSph&) const;
		int tmax;// Max number of terms retained in residual H
	public:
		/* Creator of tori in given potential. GF deemed ok if
		 * dispersion in H < tol*freqScale*Jtotal */
		TorusGenerator(const potential::BasePotential& _pot, const double _tol = 1e-9);
		Torus fitTorus(const Actions&, const double tighten = 1) const;
		Torus fitBaseTorus(const Actions&, const double tighten = 1) const;
		Torus fitFrom(const Actions&, const Torus&, const double tighten = 1) const;
		eTorus fiteTorus(const Actions&, const potential::BasePotential* _addPhi = NULL);
		eTorus fiteTorus(const Actions&, const double, const potential::BasePotential* _addPhi = NULL);
		double getDelta(double&, double&);
		double getDelta(Actions&);
		double getRsh(Actions&);
		void test_it(const Actions&, std::vector<double>&);
		std::vector<Torus> constE(const double Jrmin, const Actions& Jstart, const int Nstep);
		void old_getHn(const Torus&, int);
	};
	class EXP ActionFinderTG : public BaseActionFinder {
	private:
		const potential::PtrPotential pot;
		const ActionFinderAxisymFudge& AF;
		const TorusGenerator& TG;

	public:
		ActionFinderTG(const potential::PtrPotential& _pot, const ActionFinderAxisymFudge& _AF,
			const TorusGenerator& _TG) :
			pot(_pot), AF(_AF), TG(_TG) {
		}
		//virtual Actions actions(const coord::PosMomCyl& point) const;
		virtual ActionAngles actionAngles(const coord::PosVelCyl& point,
			Frequencies* freq = NULL) const;
		virtual ActionAngles actionAngles2(const coord::PosVelCyl& point,
			Frequencies* freq = NULL) const;
		virtual ActionAngles actionAngles3(const coord::PosVelCyl& point,
			Frequencies* freq = NULL) const;
		virtual Actions actions(const coord::PosVelCyl& point) const {
			return Actions(actionAngles(point));
		}
	};

}//namespace actions
