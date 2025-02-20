/** \file    df_halo.h
    \brief   Distribution functions for the spheroidal component (halo)
    \author  Eugene Vasiliev, James Binney
    \date    2015-2019
*/
#pragma once
#include "df_base.h"
#include "units.h"
#include "math_core.h"
#include "math_ode.h"
#include "potential_utils.h"
#include <fstream>
#include <cmath>
#define EXP __declspec(dllexport)

namespace df{

/// \name   Classes for action-based double power-law distribution function (DF)
///@{

/// Parameters that describe a double power law distribution function.
struct DoublePowerLawParam{
	double
			norm,      ///< normalization factor with the dimension of mass
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			steepness, ///< steepness of the transition between two asymptotic regimes (eta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefJrIn,  ///< contribution of radial   action to h(J), controlling anisotropy below J_0 (h_r)
			coefJzIn,  ///< contribution of vertical action to h(J), controlling anisotropy below J_0 (h_z)
			coefJrOut, ///< contribution of radial   action to g(J), controlling anisotropy above J_0 (g_r)
			coefJzOut, ///< contribution of vertical action to g(J), controlling anisotropy above J_0 (g_z)
			rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
			Jphi0,     ///< controls the steepness of rotation and the size of non-rotating core
			Jcore;     ///< central core size for a Cole&Binney-type modified double-power-law halo
	std::string Fname;
	DoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
	    coefJrIn(1), coefJzIn(1), coefJrOut(1), coefJzOut(1), rotFrac(0), Jphi0(0), Jcore(0) {}
};

/** General double power-law model.
    The distribution function is given by
    \f$  f(J) = norm / (2\pi J_0)^3
         (1 + (J_0 /h(J))^\eta )^{\Gamma / \eta}
         (1 + (g(J)/ J_0)^\eta )^{-B / \eta }
         \exp[ - (g(J) / J_{cutoff})^\zeta ]
         ( [J_{core} / h(J)]^2 - \beta J_{core} / h(J) + 1)^{\Gamma/2}   \f$,  where
    \f$  g(J) = g_r J_r + g_z J_z + g_\phi |J_\phi|  \f$,
    \f$  h(J) = h_r J_r + h_z J_z + h_\phi |J_\phi|  \f$.
    Gamma is the power-law slope of DF at small J (slopeIn), and Beta -- at large J (slopeOut),
    the transition occurs around J=J0, and its steepness is adjusted by the parameter eta.
    h_r, h_z and h_phi control the anisotropy of the DF at small J (their sum is always taken
    to be equal to 3, so that there are two free parameters -- coefJrIn = h_r, coefJzIn = h_z),
    and g_r, g_z, g_phi do the same for large J (coefJrOut = g_r, coefJzOut = g_z).
    Jcutoff is the threshold of an optional exponential suppression, and zeta measures its strength.
    Jcore is the size of the central core (if nonzero, f(J) tends to a constant limit as J --> 0
    even when the power-law slope Gamma is positive), and the auxiliary coefficient beta
    is assigned automatically from the requirement that the introduction of the core (almost)
    doesn't change the overall normalization (eq.5 in Cole&Binney 2017).
*/
class EXP DoublePowerLaw: public BaseDistributionFunction{
	const DoublePowerLawParam par;  ///< parameters of DF
	const double beta;              ///< auxiliary coefficient for the case of a central core
	public:
    /** Create an instance of double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		DoublePowerLaw(const DoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
};

/// Parameters that describe a modified double power law distribution function.
struct ModifiedDoublePowerLawParam{
	double
			norm,      ///< normalization factor with the dimension of mass
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			Jphi0,     ///< controls the steepness of rotation and the size of non-rotating core
			Jcore,     ///< central core size for a Cole&Binney-type modified double-power-law halo
			L0,        ///< helps define E surrogate jt
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefJrIn,     ///< radial anisotropy
			coefJzIn,      ///< introduces flattening
			rotFrac;   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
	std::string Fname;
	ModifiedDoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), cutoffStrength(2),
	    coefJrIn(NAN), coefJzIn(NAN), Jphi0(0), rotFrac(0), Jcore(0) {}
};

/** Modified double power-law model.
    The DF is similar to the basic double power law DF but we replace L
    by cL=a*Jz+b*Jz*|Jphi|/L+|Jphi| with a=(k+1)/2 b=(k-1)/2. When
    k=nu/Omega, the ratio of epicycle frequencies, g=Jr+.5*(1+c*xi)*cL
    is approximately a function of H, where .5*(1c*xi) is a fit to
    Omega_phi/Omegra_r
*/
class EXP ModifiedDoublePowerLaw: public BaseDistributionFunction{
	const ModifiedDoublePowerLawParam par;  ///< parameters of DF
	const double beta;              ///< auxiliary coefficient for the case of a central core
	public:
    /** Create an instance of double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		ModifiedDoublePowerLaw(const ModifiedDoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
};

/// Parameters that describe a sin double power law distribution function.
struct SinDoublePowerLawParam{
	double
			norm,      ///< normalization factor with the dimension of mass
			mass,      ///< mass
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			Jphi0,     ///< controls the steepness of rotation and the size of non-rotating core
			Jcore,     ///< central core size for a Cole&Binney-type modified double-power-law halo
			L0,        ///< helps define E surrogate jt
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			alpha,     ///< helps define xi which goes 0 -> 1 inside -> out
			beta,      ///< induces radial bias via sin
			Fin,       ///< sets coeffs a,b that determin cL
			Fout,      ///< reduces cost of low incl orbits
			rotFrac;   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
	std::string Fname;
	SinDoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(NAN), mass(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), cutoffStrength(2),
	    alpha(0.6), beta(NAN), Fin(NAN), Fout(NAN), Jphi0(0), rotFrac(0), Jcore(0) {}
};

class EXP SinDoublePowerLaw: public BaseDistributionFunction{
	const SinDoublePowerLawParam par;  ///< parameters of DF
	const double beta;              ///< auxiliary coefficient for the case of a central core
	public:
    /** Create an instance of double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		SinDoublePowerLaw(const SinDoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J) const;
		virtual void write_params(std::ofstream &strm,const units::InternalUnits &intUnits) const;
};

struct PlummerParam{
	double mass;
	double scaleRadius;
	double scaleAction;
	double mu, nu;
	PlummerParam(): mass(1), scaleRadius(1), scaleAction(1), mu(0), nu(0){}
};
class EXP PlummerDF : public BaseDistributionFunction{
	const PlummerParam par;
	double norm,Etop,Ebot,Jrtop,cLtop,Jrbot,cLbot;
	math::LinearInterpolator jrmax;
	math::LinearInterpolator cls;
	public:
		PlummerDF(const PlummerParam&);
		virtual double value(const actions::Actions& J) const;
		virtual void write_params(std::ofstream&, const units::InternalUnits&) const;
};

struct IsochroneParam{
	double mass;
	double scaleRadius;
	double mu, nu;
	IsochroneParam() : mass(1), scaleRadius(1), mu(0), nu(0) {}
};
class EXP IsochroneDF : public BaseDistributionFunction{
	const IsochroneParam par;
	double norm;
	public:
		IsochroneDF(const IsochroneParam& params);
		virtual double value(const actions::Actions& J) const;
		virtual void write_params(std::ofstream& stream,const units::InternalUnits& intUnits) const;
};
#include "Oxford.h"

///@}
}  // namespace df
