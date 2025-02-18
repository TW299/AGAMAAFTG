/** \file    actions_genfnc.h
	\brief   Generating functions for Torus mapping
	\author  Eugene Vasiliev
	\date    Feb 2016

*/
#pragma once
#include "actions_base.h"
#include "math_linalg.h"
#include <vector>
#include <cassert>
#include <map>

namespace actions {

	/** Triple index of a single term in the generating function */
	struct GenFncIndex {
		int mr, mz, mphi;
		GenFncIndex() : mr(0), mz(0), mphi(0) {}
		GenFncIndex(int _mr, int _mz, int _mphi) : mr(_mr), mz(_mz), mphi(_mphi) {}
	};

	/** Array of indices that represent all non-trivial terms in the generating function */
	typedef std::vector<GenFncIndex> GenFncIndices;
	typedef std::vector<Actions> GenFncDerivs;
	struct GenFncFrac {
		int mz, mphi, krmin;
		GenFncFrac() : mz(0), mphi(0) {}
		GenFncFrac(int _mz, int _mphi = 0, int _krmin = 0) : mz(_mz), mphi(_mphi), krmin(_krmin) {}
	};

	typedef std::vector<GenFncFrac> GenFncFracs;
	/** Generating function that maps the true action/angles to toy action/angles */
	class EXP GenFnc {
	public:
		GenFnc() {}
		GenFnc(const GenFncIndices _indices, const std::vector<double> _values,
			const GenFncDerivs _derivs) : indices(_indices), values(_values), derivs(_derivs) {
			
		}
		GenFnc(const GenFncIndices _indices, const std::vector<double> _values,
			const GenFncDerivs _derivs,const GenFncFracs _fracs) : indices(_indices), values(_values), derivs(_derivs),fracs(_fracs){
		}
		GenFnc(const GenFnc& GF) : indices(GF.indices), values(GF.values), derivs(GF.derivs),fracs(GF.fracs) {}
		void print(void) const;
		void write(FILE*) const;
		void read(FILE*);
		double giveValue(const GenFncIndex&) const;
		unsigned int numParams() const { return indices.size(); }
		std::pair<int, int> maxIndices() const;
		Actions toyJ(const Actions&, const Angles&) const;// returns toy actions
		Angles trueA(const Angles& thetaT) const;// returns true angles
		Angles  toyA(const Angles& theta) const;//returns toy angles
		ActionAngles true2toy(const ActionAngles& aa) const;//input true aa, returns aaT
		DerivAng<coord::Cyl> dJdt(const Angles&) const;//d JT/dthetaT at fixed J
		double dtbydtT_Jacobian(const Angles&, math::Matrix<double>&) const;
		double TermRatio(int kstart, int kz, double& B,double *r=NULL) const;
		GenFnc reSum(int kz, const double Smin, double& B, double& b) const;
		GenFnc operator * (const double a) const;
		GenFnc& operator *= (const double a);
		GenFnc& operator += (const GenFnc& GF);
		GenFnc operator + (const GenFnc& GF) const;
		//	private:
		GenFncIndices indices;      ///< indices of terms in generating function
		std::vector<double> values; ///< amplitudes of terms in generating function
		GenFncDerivs derivs;        ///< amplitudes of derivatives dS/dJ_{r,z,phi}
		GenFncFracs fracs;
	};

	/** Variant of generating function used during the fitting process.
		It converts true actions to toy actions at any of the pre-determined array of angles
		specified at its construction, for the given amplitudes of its terms;
		and also provides the derivatives of toy actions by these amplitudes.
		Unlike the general-purpose action map, it only operates on the restricted set of angles,
		which avoids repeated computation of trigonometric functions; on the other hand,
		the amplitudes of its terms are not fixed at construction, but provided at each call.
	*/
	class GenFncFit {
	public:
		/** construct the object for the fixed indexing scheme,
			values of real actions, and array of toy angles */
		GenFncFit(const GenFncIndices& indices, const Actions& acts);
		GenFncFit(void) {};
		/** number of terms in the generating function */
		unsigned int numParams() const { return indices.size(); }

		/** number of points in the array of angles */
		unsigned int numPoints() const { return angs.size(); }

		/** perform mapping from real actions to toy actions at the specific values of angles.
			\param[in]  indexAngle is the index of element in the pre-determined grid of angles;
			\param[in]  values     is the array of amplitudes of the generating function terms;
			\return  toy actions and angles.
		*/
		ActionAngles toyActionAngles(const unsigned int indexAngle, const double values[]) const;

		GenFncIndices expand(std::vector<double>& params);
		/** compute derivatives of toy actions w.r.t generating function coefficients.
			\param[in]  indexAngle is the index of element in the pre-determined grid of angles;
			\param[in]  indexCoef  is the index of term in the generating function;
			\return  derivatives of each action by the amplitude of this term.
		*/
		inline Actions deriv(const unsigned int indexAngle, const unsigned int indexCoef) const {
			double val = coefs(indexAngle, indexCoef);  // no range check performed!
			return Actions(
				val * indices[indexCoef].mr,
				val * indices[indexCoef].mz,
				val * indices[indexCoef].mphi);
		}
		void reset(const Actions J) {
			acts = J;
		}
		void print(const std::vector<double>& params) const;
		void write(FILE*, const std::vector<double>& params) const;
	private:
		GenFncIndices indices;    ///< indices of terms in generating function
		Actions acts;             ///< values of real actions
		std::vector<Angles> angs; ///< grid of toy angles
		math::Matrix<double> coefs;     ///< precomputed trigonometric functions at the grid of angles
	};

	


	class EXP GenFncFitSeries {
	public:
		/** construct the object for the fixed indexing scheme,
			values of real actions, and array of toy angles */
		GenFncFitSeries(const GenFncIndices& indices, const GenFncFracs& fracs,
			const Actions& acts);
		GenFncFitSeries(void) {};
		/** number of terms in the generating function */
		unsigned int numTerms() const { return indices.size(); }

		/** number of fractions in the generating function */
		unsigned int numFracs() const { return fracs.size(); }

		/** number of parameters in the generating function */
		unsigned int numParams() const { return indices.size() + 2 * fracs.size(); }

		/** number of points in the array of angles */
		unsigned int numPoints() const { return angs.size(); }

		std::pair<int, int> maxIndices() const;

		/** perform mapping from real actions to toy actions at the specific values of angles.
		\param[in]  indexAngle is the index of element in the pre-determined grid of angles;
		\param[in]  values     is the array of amplitudes of the generating function terms;
		\return  toy actions and angles.
	*/
		ActionAngles toyActionAngles(const unsigned int indexAngle, const double values[]) const;

		GenFncIndices expand(std::vector<double>& params);

		/** compute derivatives of toy actions w.r.t generating function coefficients.
		 \param[in]  indexAngle is the index of element in the pre-determined grid of angles;
		 \param[in]  indexCoef  is the index of parameter in the generating function;
		 \return  derivatives of each action by the parameter.
	 */
		Actions deriv(const unsigned int indexAngle, const unsigned int indexCoef, const double[]) const;

		void reset(const Actions J) {
			acts = J;
		}
		void print(const std::vector<double>& params) const;
		void write(FILE*, const std::vector<double>& params) const;
		void read(FILE*, std::vector<double>& params);
		GenFncIndices indices;    ///< indices of terms in generating function
		GenFncFracs fracs;        ///< fractions summing series in theta_r
		Actions acts;             ///< values of real actions
	private:
		std::vector<Angles> angs; ///< grid of toy angles
		math::Matrix<double> coefs;     ///< precomputed trigonometric functions at the grid of angles
	};


}  // namespace actions
