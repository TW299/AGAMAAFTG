#include "actions_newgenfnc.h"
#include "math_fit.h"
#include "math_core.h"
#include <cmath>

namespace actions {

	namespace { // internal

		/// create an array of angles uniformly covering the range [0:pi]  (NB: why not 2pi?)
		static std::vector<Angles> makeGridAngles(unsigned int nr, unsigned int nz, unsigned int nphi = 1)
		{
			std::vector<Angles> vec(nr * nz * nphi);
			for (unsigned int ir = 0; ir < nr; ir++) {
				double thetar = ir * M_PI / nr;
				for (unsigned int iz = 0; iz < nz; iz++) {
					double thetaz = iz * M_PI / nz;
					for (unsigned int iphi = 0; iphi < nphi; iphi++)
						vec[(ir * nz + iz) * nphi + iphi] =
						Angles(thetar, thetaz, iphi * M_PI / fmax(nphi - 1, 1));
				}
			}
			//printf("Grid of %d angles created\n", nr * nz * nphi);
			return vec;
		}

		/// create grid in angles with size determined by the maximal Fourier harmonic in the indices array
		static std::vector<Angles> makeGridAngles(const GenFncIndices& indices)
		{
			int maxmr = 4, maxmz = 4, maxmphi = 0;
			for (unsigned int i = 0; i < indices.size(); i++) {
				maxmr = std::max<int>(maxmr, math::abs(indices[i].mr));
				maxmz = std::max<int>(maxmz, math::abs(indices[i].mz));
				maxmphi = std::max<int>(maxmphi, math::abs(indices[i].mphi));
			}
			//maxmz /= 2;
			return makeGridAngles(6 * (maxmr / 4 + 1), 6 * (maxmz / 4 + 1), maxmphi > 0 ? 6 * (maxmphi / 4 + 1) : 1);
		}

		/// return the absolute value of an element in a map, or zero if it doesn't exist
		static inline double absvalue(const std::map< std::pair<int, int>, double >& indPairs, int ir, int iz)
		{
			if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
				return fabs(indPairs.find(std::make_pair(ir, iz))->second);
			else
				return 0;
		}

		/** Helper class to be used in the iterative solution of a nonlinear system of equations
			that implicitly define the toy angles as functions of real angles and the derivatives
			of the generating function.
		*/
		class AngleFinder : public math::IFunctionNdimDeriv {
		public:
			/*AngleFinder(const GenFncIndices& _indices, const GenFncDerivs& _dSby, const Angles& _ang) :
				indices(_indices), dSby(_dSby), ang(_ang) {
				
			}*/
			AngleFinder(const GenFncIndices& _indices, const GenFncDerivs& _dSby, const Angles& _ang, const std::vector<double> _params, const GenFncFracs _fracs) :
				indices(_indices), dSby(_dSby), ang(_ang), params(_params), fracs(_fracs) {
			}
			virtual unsigned int numVars() const { return 3; }
			virtual unsigned int numValues() const { return 3; }

			virtual void evalDeriv(const double vars[], double values[], double* derivs = 0) const
			{
				if (values) {
					values[0] = vars[0] - ang.thetar;
					values[1] = vars[1] - ang.thetaz;
					values[2] = vars[2] - ang.thetaphi;
				}
				if (derivs) {
					derivs[0] = derivs[4] = derivs[8] = 1.;  // diagonal
					derivs[1] = derivs[2] = derivs[3] = derivs[5] = derivs[6] = derivs[7] = 0;  // off-diag
				}
				for (unsigned int i = 0; i < indices.size(); i++) {
					if (i < indices.size()) {
						double arg = indices[i].mr * vars[0] + indices[i].mz * vars[1] +
							indices[i].mphi * vars[2];    // argument of trig functions
						if (values) {
							double s = sin(arg);
							values[0] += s * dSby[i].Jr;
							values[1] += s * dSby[i].Jz;
							values[2] += s * dSby[i].Jphi;
						}
						if (derivs) {
							double c = cos(arg);
							derivs[0] += c * dSby[i].Jr * indices[i].mr;
							derivs[1] += c * dSby[i].Jr * indices[i].mz;
							derivs[2] += c * dSby[i].Jr * indices[i].mphi;
							derivs[3] += c * dSby[i].Jz * indices[i].mr;
							derivs[4] += c * dSby[i].Jz * indices[i].mz;
							derivs[5] += c * dSby[i].Jz * indices[i].mphi;
							derivs[6] += c * dSby[i].Jphi * indices[i].mr;
							derivs[7] += c * dSby[i].Jphi * indices[i].mz;
							derivs[8] += c * dSby[i].Jphi * indices[i].mphi;
						}
					}
				}
				unsigned int fracNum = 0;
				for (unsigned int indexCoef = indices.size(); indexCoef < params.size(); indexCoef += 2) {
					Actions dBdJ = dSby[indexCoef], dpsidJ = dSby[indexCoef + 1];
					int kmin = fracs[fracNum].krmin;
					const int mz = fracs[fracNum].mz;
					double B = params[indexCoef];
					double b = tanh(params[indexCoef + 1]);
					std::complex<double> tz(exp(std::complex<double>
						(0, mz * vars[1])));
					std::complex<double> tr(exp(std::complex<double>
						(0, vars[0])));
					std::complex<double> trk(exp(std::complex<double>
						(0, kmin * vars[0])));
					std::complex<double> onembtr((std::complex<double>(1, 0) - b * tr));

					double F = imag(tz * trk / onembtr);
					double dFdthetar = real(b * tr * trk * tz / pow_2(onembtr));
					double dFdthetaz = real(double(mz) * trk * tz / (onembtr));
					double dFdthetarphi = mz * (1 - b * b) * real(trk * tz / pow_2(onembtr));
					double dFdthetazphi = (1 - b * b) * real(trk * tz / pow_2(onembtr) * (1.0 + 2 * b * tr / onembtr + double(kmin)));
					double dFdphi = (1.0 - b * b) * imag(tz * tr * trk / pow_2(onembtr));
					if (values) {
						values[0] += B * dFdphi * dpsidJ.Jr + F * dBdJ.Jr;
						values[1] += B * dFdphi * dpsidJ.Jz + F * dBdJ.Jz;
						values[2] += B * dFdphi * dpsidJ.Jphi + F * dBdJ.Jphi;
					}if (derivs) {
						derivs[0] += B * dFdthetarphi * dpsidJ.Jr + dFdthetar * dBdJ.Jr;
						derivs[1] += B * dFdthetazphi * dpsidJ.Jr + dFdthetaz * dBdJ.Jr;

						derivs[3] += B * dFdthetarphi * dpsidJ.Jz + dFdthetar * dBdJ.Jz;
						derivs[4] += B * dFdthetazphi * dpsidJ.Jz + dFdthetaz * dBdJ.Jz;

						derivs[6] += B * dFdthetarphi * dpsidJ.Jphi + dFdthetar * dBdJ.Jphi;
						derivs[7] += B * dFdthetazphi * dpsidJ.Jphi + dFdthetaz * dBdJ.Jphi;
					}
					fracNum++;
				}
			}
		private:
			const GenFncIndices& indices; ///< indices of terms in generating function
			const GenFncDerivs& dSby;     ///< amplitudes of derivatives dS/dJ_{r,z,phi}
			const GenFncFracs& fracs;
			const std::vector<double>& params;
			const Angles ang;             ///< true angles
		};

	} // internal namespace

	void GenFnc::write(FILE* ofile) const {
		fprintf(ofile, "%zd\n", indices.size());
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g %g %g %g\n",
				indices[i].mr, indices[i].mz, indices[i].mphi,
				values[i], derivs[i].Jr, derivs[i].Jz,
				derivs[i].Jphi);
	}
	void GenFnc::read(FILE* ifile) {
		int n;
		fscanf_s(ifile, "%d", &n);
		indices.resize(n); values.resize(n); derivs.resize(n);
		for (int i = 0; i < indices.size(); i++)
			fscanf_s(ifile, "%d %d %d %g %g %g %g\n",
				&indices[i].mr, &indices[i].mz, &indices[i].mphi,
				&values[i], &derivs[i].Jr, &derivs[i].Jz, &derivs[i].Jphi);
	}
	void GenFnc::print(void) const {
		for (int i = 0; i < values.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, values[i]);
		printf("\n");
	}

	Actions GenFnc::toyJ(const Actions& J, const Angles& thetaT) const {//returns toyJ(thetaT)
		Actions JT(J);
		for (unsigned int i = 0; i < indices.size(); i++) {// compute toy actions
			double val = values[i] * cos(indices[i].mr * thetaT.thetar +
				indices[i].mz * thetaT.thetaz +
				indices[i].mphi * thetaT.thetaphi);
			JT.Jr += val * indices[i].mr;
			JT.Jz += val * indices[i].mz;
			JT.Jphi += val * indices[i].mphi;
		}unsigned int fracNum = 0;
		for (unsigned int j = indices.size();j < values.size();j+=2) {
			double B = values[j], b = tanh(values[j + 1]);
			int kmin = fracs[fracNum].krmin;
			const int mz = fracs[fracNum].mz;
			std::complex<double> tz(exp(std::complex<double>
				(0, mz * thetaT.thetaz)));
			std::complex<double> tr(exp(std::complex<double>
				(0, thetaT.thetar)));
			std::complex<double> trk(exp(std::complex<double>
				(0, kmin * thetaT.thetar)));
			std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
			JT.Jr += B * real(trk * tz * (b * tr / pow_2(onembtr) + double(kmin) / onembtr));
			JT.Jz += B * mz * real(tz * trk / onembtr);
			fracNum++;
		}
		if (JT.Jr < 0 || JT.Jz < 0) {
			//printf("Ji<0 %f %f %f %f\n",JT.Jr,JT.Jz, thetaT.thetar, thetaT.thetaz);
			JT.Jr = fmax(JT.Jr, 0); JT.Jz = fmax(JT.Jz, 0);    // prevent non-physical negative values
		}
		return JT;
	}
	Angles GenFnc::trueA(const Angles& thetaT) const {// toy Angs -> real Angs
		Angles theta(thetaT);
		for (unsigned int i = 0; i < indices.size(); i++) {
			double sinT = sin(indices[i].mr * thetaT.thetar + indices[i].mz * thetaT.thetaz
				+ indices[i].mphi * thetaT.thetaphi);
			theta.thetar += derivs[i].Jr * sinT;
			theta.thetaz += derivs[i].Jz * sinT;
			theta.thetaphi += derivs[i].Jphi * sinT;
		}if (fracs.size() > 0) {
			for (unsigned int i = indices.size(); i < fracs.size(); i++) {
				Actions dBdJ = derivs[i], dpsidJ = derivs[i + 1];
				int kmin = fracs[i - indices.size()].krmin;
				const int mz = fracs[i - indices.size()].mz;
				double B = values[i];double b = tanh(values[i + 1]);
				
				std::complex<double> tz(exp(std::complex<double>
					(0, mz * thetaT.thetaz)));
				std::complex<double> tr(exp(std::complex<double>
					(0, thetaT.thetar)));
				std::complex<double> trk(exp(std::complex<double>
					(0, kmin * thetaT.thetar)));
				std::complex<double> onembtr((std::complex<double>(1, 0) - b * tr));

				double F = imag(tz * trk / onembtr);
				double dFdthetar = real(b * tr * trk * tz / pow_2(onembtr));
				double dFdthetaz = real(double(mz) * trk * tz / (onembtr));
				double dFdthetarphi = mz * (1 - b * b) * real(trk * tz / pow_2(onembtr));
				double dFdthetazphi = (1 - b * b) * real(trk * tz / pow_2(onembtr) * (1.0 + 2 * b * tr / onembtr + double(kmin)));
				double dFdphi = (1.0 - b * b) * imag(tz * tr * trk / pow_2(onembtr));
				theta.thetar+= B * dFdphi * dpsidJ.Jr + F * dBdJ.Jr;
				theta.thetaz += B * dFdphi * dpsidJ.Jz + F * dBdJ.Jz;
				theta.thetaz+= B * dFdphi * dpsidJ.Jphi + F * dBdJ.Jphi;
			}
		}
		return theta;
	}

	Angles GenFnc::toyA(const Angles& theta) const {//real Angs -> toy Angs
		double thetaT[3], Theta[3] = { theta.thetar, theta.thetaz, theta.thetaphi };
		AngleFinder AF(indices, derivs, theta,values,fracs);
		int numIter = math::findRootNdimDeriv(AF, Theta, 1e-6, 10, thetaT);
		if (numIter >= 10) printf("Max iterations in findRootNdimDeriv\n");
		return Angles(math::wrapAngle(thetaT[0]),
			math::wrapAngle(thetaT[1]), math::wrapAngle(thetaT[2]));
	}

	ActionAngles GenFnc::true2toy(const ActionAngles& aa) const {
		Angles thetaT(toyA(Angles(aa)));
		Actions JT(toyJ(Actions(aa), thetaT));
		return ActionAngles(JT, thetaT);
	}

	DerivAng<coord::Cyl> GenFnc::dJdt(const Angles& thetaT) const {
		DerivAng<coord::Cyl> D;
		for (unsigned int i = 0; i < indices.size(); i++) {
			D.dbythetar.R = 0; D.dbythetaz.R = 0; D.dbythetaphi.R = 0;
			D.dbythetar.z = 0; D.dbythetaz.z = 0; D.dbythetaphi.z = 0;
			D.dbythetar.phi = 0; D.dbythetaz.phi = 0; D.dbythetaphi.phi = 0;
		}
		for (unsigned int i = 0; i < indices.size(); i++) {
			double sinT = values[i] * sin(indices[i].mr * thetaT.thetar
				+ indices[i].mz * thetaT.thetaz
				+ indices[i].mphi * thetaT.thetaphi);
			D.dbythetar.R -= indices[i].mr * indices[i].mr * sinT;//derivs of Jr
			D.dbythetaz.R -= indices[i].mr * indices[i].mz * sinT;
			D.dbythetaphi.R -= indices[i].mr * indices[i].mphi * sinT;
			D.dbythetar.z -= indices[i].mz * indices[i].mr * sinT;//derivs of Jz
			D.dbythetaz.z -= indices[i].mz * indices[i].mz * sinT;
			D.dbythetaphi.z -= indices[i].mz * indices[i].mphi * sinT;
			D.dbythetar.phi -= indices[i].mphi * indices[i].mr * sinT;//derivs of Jphi
			D.dbythetaz.phi -= indices[i].mphi * indices[i].mz * sinT;
			D.dbythetaphi.phi -= indices[i].mphi * indices[i].mphi * sinT;
		}
		return D;
	}

	double GenFnc::dtbydtT_Jacobian(const Angles& thetaT, math::Matrix<double>& M) const {// dtheta_i/dthetaT_j
		M(0, 0) = M(1, 1) = M(2, 2) = 1; M(0, 1) = M(0, 2) = M(1, 2) = M(1, 0) = M(2, 0) = M(2, 1) = 0;
		for (unsigned int i = 0; i < indices.size(); i++) {
			double cosT = cos(indices[i].mr * thetaT.thetar
				+ indices[i].mz * thetaT.thetaz
				+ indices[i].mphi * thetaT.thetaphi);
			M(0, 0) += derivs[i].Jr * indices[i].mr * cosT;
			M(0, 1) += derivs[i].Jr * indices[i].mz * cosT;
			M(0, 2) += derivs[i].Jr * indices[i].mphi * cosT;
			M(1, 0) += derivs[i].Jz * indices[i].mr * cosT;
			M(1, 1) += derivs[i].Jz * indices[i].mz * cosT;
			M(1, 2) += derivs[i].Jz * indices[i].mphi * cosT;
			M(2, 0) += derivs[i].Jphi * indices[i].mr * cosT;
			M(2, 1) += derivs[i].Jphi * indices[i].mz * cosT;
			M(2, 2) += derivs[i].Jphi * indices[i].mphi * cosT;
		}
		double det = M(0, 0) * (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2))
			- M(1, 0) * (M(0, 1) * M(2, 2) - M(2, 1) * M(0, 2))
			+ M(2, 0) * (M(0, 1) * M(1, 2) - M(1, 1) * M(0, 2));
		return fabs(det);
	}

	double GenFnc::giveValue(const GenFncIndex& ind) const {
		for (int i = 0; i < indices.size(); i++) {
			if (ind.mr != indices[i].mr) continue;
			if (ind.mz != indices[i].mz) continue;
			if (ind.mphi == indices[i].mphi) return values[i];
		}
		return NAN;
	}

	double GenFnc::TermRatio(int kStart, int kz, double& B,double *dispersion) const {
		double Slast, rbar = 0; int nbar = 0;
		for (int i = kStart; i < maxIndices().first; i++) {//get term ratio
			double S = giveValue(GenFncIndex(i, kz, 0));
			if (std::isnan(S) || i - kStart > 20) break;
			if (i == kStart) {
				B = S;
			}
			else{
				
				rbar += (S / Slast); nbar++;
			}
			Slast = S;
		}
		double disp = 0;
		double b = rbar / nbar;
		B /= pow(b, kStart);
		if (dispersion) {
			for (int i = kStart; i < maxIndices().first; i++) {//get term ratio
				double S = giveValue(GenFncIndex(i, kz, 0));
				if (std::isnan(S) || i - kStart > 20) break;
				disp += pow_2(B * pow(b, (i)) - S) / nbar;
			}
			*dispersion = sqrt(disp);
		}
		return b;
	}

	GenFnc GenFnc::reSum(int kz, const double Smin, double& B, double& b) const {
		std::pair<int, int> maxIs(maxIndices());
		int kstart = 4;
		int k1 = kstart;
		b = TermRatio(kstart, kz, B);
		GenFncIndices newIndices;
		std::vector<double> newValues;
		GenFncDerivs newDerivs;
		for (int i = 0; i < indices.size(); i++) {
			if (indices[i].mz == kz&& indices[i].mr >= k1) {
				double S = values[i] - B * pow(b, indices[i].mr);
				if (fabs(S) > Smin) {
					newIndices.push_back(indices[i]);
					newValues.push_back(S);
					if (derivs.size() > 0) {
						newDerivs.push_back(derivs[i]);//subtraction required!!
					}
				}
			}
			else {
				newIndices.push_back(indices[i]);
				newValues.push_back(values[i]);
				if (derivs.size() > 0) {
					newDerivs.push_back(derivs[i]);//subtraction required!!
				}
			}
		}
		return GenFnc(newIndices, newValues, newDerivs);
	}

	std::pair<int, int> GenFnc::maxIndices() const {
		int mr = 0, mz = 0;
		for (int i = 0; i < indices.size(); i++) {
			mr = std::max<int>(mr, std::abs(indices[i].mr));
			mz = std::max<int>(mz, std::abs(indices[i].mz));
		}
		return std::make_pair(mr, mz);
	}


	GenFnc GenFnc::operator * (const double a) const {
		std::vector<double> values2(values.size());
		for (int i = 0; i < values.size(); i++)
			values2[i] = a * values[i];
		std::vector<Actions> derivs2(derivs.size());
		for (int i = 0; i < derivs.size(); i++)
			derivs2[i] = derivs[i] * a;
		GenFnc G2(indices, values2, derivs2);
		return G2;
	}
	GenFnc& GenFnc::operator *= (const double a) {
		for (int i = 0; i < values.size(); i++)
			values[i] *= a;
		for (int i = 0; i < derivs.size(); i++)
			derivs[i] *= a;
		return *this;
	}
	GenFnc& GenFnc::operator += (const GenFnc& G) {
		GenFnc G2(G);
		for (int i = 0; i < indices.size(); i++) {
			int mr = indices[i].mr, mz = indices[i].mz, mphi = indices[i].mphi;
			std::vector<double>::iterator jt = G2.values.begin();
			std::vector<Actions>::iterator kt = G2.derivs.begin();//find terms matching present ones
			for (GenFncIndices::iterator it = G2.indices.begin(); it != G2.indices.end();) {
				if (mr == (*it).mr && mz == (*it).mz && mphi == (*it).mphi) {
					values[i] += (*jt); derivs[i] += (*kt);
					it = G2.indices.erase(it);
					jt = G2.values.erase(jt);
					kt = G2.derivs.erase(kt);
					break;
				}
				else {
					it++; jt++; kt++;
				}
			}
			//printf("(%zd %zd %zd) ",G2.indices.size(),G2.values.size(),G2.derivs.size());
		}
		if (G2.indices.size() > 0) {
			for (int j = 0; j < G2.indices.size(); j++) {
				indices.push_back(G2.indices[j]);
				values.push_back(G2.values[j]);
				derivs.push_back(G2.derivs[j]);
			}
		}
		return *this;
	}

	GenFnc GenFnc::operator + (const GenFnc& GF) const {
		GenFnc G2 = *this;
		G2 += GF;
		return G2;
	}

	GenFncFit::GenFncFit(const GenFncIndices& _indices, const Actions& _acts) :
		indices(_indices), acts(_acts)
	{
		angs = makeGridAngles(indices);
		coefs = math::Matrix<double>(angs.size(), indices.size());

		for (unsigned int indexAngle = 0; indexAngle < angs.size(); indexAngle++)
			for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++)
				coefs(indexAngle, indexCoef) = cos(
					indices[indexCoef].mr * angs[indexAngle].thetar +
					indices[indexCoef].mz * angs[indexAngle].thetaz +
					indices[indexCoef].mphi * angs[indexAngle].thetaphi);
	}

	ActionAngles GenFncFit::toyActionAngles(unsigned int indexAngle, const double values[]) const
	{
		ActionAngles aa(acts, angs[indexAngle]);
		for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++) {
			double val = values[indexCoef] * coefs(indexAngle, indexCoef);
			aa.Jr += val * indices[indexCoef].mr;
			aa.Jz += val * indices[indexCoef].mz;
			aa.Jphi += val * indices[indexCoef].mphi;
		}
		// non-physical negative actions may appear,
		// which means that these values of parameters are unsuitable.
		return aa;
	}
	void GenFncFit::print(const std::vector<double>& params) const {
		for (int i = 0; i < params.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
		printf("\n");
	}
	void GenFncFit::write(FILE* ofile, const std::vector<double>& params) const {
		fprintf(ofile, "%zd\n", indices.size());
		for (int i = 0; i < params.size(); i++)
			fprintf(ofile, "%d %d %d %g\n", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
	}

	GenFncIndices GenFncFit::expand(std::vector<double>& params)
	{   /// NOTE: here we specialize for the case of axisymmetric systems!
		assert(params.size() == numParams());
		std::map< std::pair<int, int>, double > indPairs;
		int numadd = 0;
		GenFncIndices newIndices(indices);

		// 1. Store amplitudes & determine the extent of existing grid in (mr,mz)
		int maxmr = 0, maxmz = 0;
		for (unsigned int i = 0; i < indices.size(); i++) {
			indPairs[std::make_pair(indices[i].mr, indices[i].mz)] = params[i];
			maxmr = std::max<int>(maxmr, math::abs(indices[i].mr));
			maxmz = std::max<int>(maxmz, math::abs(indices[i].mz));
		}
		/*if(maxmz==0) {  // dealing with the case Jz==0 -- add only two elements in m_r
			newIndices.push_back(GenFncIndex(maxmr+1, 0, 0));
			newIndices.push_back(GenFncIndex(maxmr+2, 0, 0));
			return newIndices;
		}*/

		// 2. determine the largest amplitude of coefs that are at the boundary of existing values
		double maxval = 0;
		int indexz = 0;
		//note terms that are not themselves off end, but have a neighbour that is
		for (int ir = 0; ir <= maxmr + 2; ir++)
			for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end() &&
					((iz <= 0 && indPairs.find(std::make_pair(ir, iz - 2)) == indPairs.end()) ||
						(iz >= 0 && indPairs.find(std::make_pair(ir, iz + 2)) == indPairs.end()) ||
						indPairs.find(std::make_pair(ir + 1, iz)) == indPairs.end()))
					maxval = fmax(fabs(indPairs[std::make_pair(ir, iz)]), maxval);
				if( fabs(indPairs[std::make_pair(ir, iz)]) > maxval){
					indexz = iz;
				}
			}
		printf("Index %d maxval %g\n",indexz, maxval);

		// 4. add more terms adjacent to the existing ones at the boundary, if they are large enough
		double thresh = maxval * 0.1;
		for (int ir = 0; ir <= maxmr + 3; ir++)
			for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end()
					|| (ir == 0 && iz >= 0)) continue;  // already exists or not required
				if (absvalue(indPairs, ir - 3, iz) >= thresh ||
					absvalue(indPairs, ir - 2, iz) >= thresh ||
					absvalue(indPairs, ir - 1, iz) >= thresh ||
					absvalue(indPairs, ir, iz - 2) >= thresh ||
					absvalue(indPairs, ir, iz + 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz - 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz + 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz) >= thresh)
				{   // add a term if any of its neighbours are large enough
					newIndices.push_back(GenFncIndex(ir, iz, 0));
					numadd++;
				}
			}

		// 4. finish up
		printf("%d terms added overall\n", numadd);
		assert(numadd > 0);
		return newIndices;
	}

	GenFncFitSeries::GenFncFitSeries(const GenFncIndices& _indices, const GenFncFracs& _fracs, const Actions& _acts) :
		indices(_indices), fracs(_fracs), acts(_acts)
	{
		angs = makeGridAngles(indices);
		coefs = math::Matrix<double>(angs.size(), indices.size());

		for (unsigned int indexAngle = 0; indexAngle < angs.size(); indexAngle++)
			for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++)
				coefs(indexAngle, indexCoef) = cos(
					indices[indexCoef].mr * angs[indexAngle].thetar +
					indices[indexCoef].mz * angs[indexAngle].thetaz +
					indices[indexCoef].mphi * angs[indexAngle].thetaphi);
	}

	std::pair<int, int> GenFncFitSeries::maxIndices() const {
		int mr = 0, mz = 0;
		for (int i = 0; i < indices.size(); i++) {
			mr = std::max<int>(mr, std::abs(indices[i].mr));
			mz = std::max<int>(mz, std::abs(indices[i].mz));
		}
		return std::make_pair(mr, mz);
	}


	/* To handle long series in theta_r we add
	 * exp(i*m*theta_z)/(1-b*exp(i*theta_r)) with b=tanh(psi) to keep |b|<1
	*/
	ActionAngles GenFncFitSeries::toyActionAngles(unsigned int indexAngle, const double values[]) const
	{
		ActionAngles aa(acts, angs[indexAngle]);
		for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++) {
			double val = values[indexCoef] * coefs(indexAngle, indexCoef);
			aa.Jr += val * indices[indexCoef].mr;
			aa.Jz += val * indices[indexCoef].mz;
			aa.Jphi += val * indices[indexCoef].mphi;
		}
		unsigned int fracNum = 0;
		
		for (unsigned int indexCoef = indices.size(); indexCoef < numParams(); indexCoef += 2) {
			double B = values[indexCoef], b = tanh(values[indexCoef + 1]);
			int kmin = fracs[fracNum].krmin;
			const int mz = fracs[fracNum].mz;
			std::complex<double> tz(exp(std::complex<double>
				(0, mz * angs[indexAngle].thetaz)));
			std::complex<double> tr(exp(std::complex<double>
				(0, angs[indexAngle].thetar)));
			std::complex<double> trk(exp(std::complex<double>
				(0, kmin*angs[indexAngle].thetar)));
			std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
			aa.Jr += B *  real(trk*tz*(b*tr/ pow_2(onembtr)+double(kmin)/onembtr));
			aa.Jz += B * mz * real(tz*trk / onembtr);
			fracNum++;
		}
		return aa;
	}
	Actions GenFncFitSeries::deriv(const unsigned int indexAngle, const unsigned int indexCoef,
		const double values[]) const {
		Actions dJ;
		if (indexCoef < indices.size()) {//Dif wrt S_k
			double val = coefs(indexAngle, indexCoef);  // no range check performed!
			dJ = Actions(
				val * indices[indexCoef].mr,
				val * indices[indexCoef].mz,
				val * indices[indexCoef].mphi);
		}
		else {//dif wrt parameter of a fraction
			int fracNum = (indexCoef - indices.size()) / 2, k = (indexCoef - indices.size()) % 2;
			const double B = values[indices.size() + 2 * fracNum],
				b = tanh(values[indices.size() + 2 * fracNum + 1]), onembsq = 1 - b * b;
			int kmin = fracs[fracNum].krmin;
			const int mz = fracs[fracNum].mz;
			std::complex<double> tz(exp(std::complex<double>
				(0, mz * angs[indexAngle].thetaz)));
			std::complex<double> tr(exp(std::complex<double>
				(0, angs[indexAngle].thetar)));
			std::complex<double> trk(exp(std::complex<double>
				(0, kmin*angs[indexAngle].thetar)));
			std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
			if (k == 0) {//dif wrt B
				dJ.Jz = mz * real(tz*trk / onembtr);
				dJ.Jr = real(trk * tz * (b * tr / pow_2(onembtr) + double(kmin) / onembtr));

			}
			//aa.Jr += B *  real(trk*tz*(b*tr/ pow_2(onembtr)+double(kmin)/onembtr));
			//aa.Jz += B * mz * real(tz * trk / onembtr);
			else {//dif wrt b
				std::complex<double> Z(tr*tz*trk / pow_2(onembtr));
				dJ.Jz = mz * B * onembsq * real(Z);
				dJ.Jr = B * onembsq * real(Z * (1.0+ 2 * b * tr/onembtr+ double(kmin)));
			}
			dJ.Jphi = 0;
		}
		return dJ;
	}
	void GenFncFitSeries::print(const std::vector<double>& params) const {
		printf("Fourier oeffs\n");
		for (int i = 0; i < indices.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
		printf("\n");
		printf("Fractions\n");
		for (int i = 0; i < fracs.size(); i++)
			printf("m: %d, B: %g, b: %g\n",
				fracs[i].mz, params[indices.size() + 2 * i], params[indices.size() + 2 * i + 1]);
	}
	void GenFncFitSeries::write(FILE* ofile, const std::vector<double>& params) const {
		fprintf(ofile, "%zd\n", indices.size());
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g\n", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
		fprintf(ofile, "%zd\n", fracs.size());
		for (int i = 0; i < fracs.size(); i++)
			fprintf(ofile, "%d %g %g\n",//mz, B, b
				fracs[i].mz, params[indices.size() + 2 * i], params[indices.size() + 2 * i + 1]);
	}
	void GenFncFitSeries::read(FILE* ifile, std::vector<double>& params) {
		indices.clear(); params.clear();
		int mr, mz, mphi; double p;
		int sI; fscanf_s(ifile, "%d", &sI);
		for (int i = 0; i < sI; i++) {
			fscanf_s(ifile, "%d %d %d %lg", &mr, &mz, &mphi, &p);
			indices.push_back(GenFncIndex(mr, mz, mphi));
			params.push_back(p);
		}
		printf("Read %zd terms\n", params.size());
		int nf; double B, b;
		fscanf_s(ifile, "%d", &nf);
		for (int i = 0; i < nf; i++) {
			fscanf_s(ifile, "%d %lg %lg", &mz, &B, &b);
			fracs.push_back(GenFncFrac(mz));
			params.push_back(B); params.push_back(b);
		}
		printf("Read %zd fracs\n", fracs.size());
	}

	GenFncIndices GenFncFitSeries::expand(std::vector<double>& params)
	{   /// NOTE: here we specialize for the case of axisymmetric systems!
		assert(params.size() == numParams());
		std::map< std::pair<int, int>, double > indPairs;
		int numadd = 0;
		GenFncIndices newIndices(indices);

		// 1. Store amplitudes & determine the extent of existing grid in (mr,mz)
		int maxmr = 0, maxmz = 0;
		for (unsigned int i = 0; i < indices.size(); i++) {
			indPairs[std::make_pair(indices[i].mr, indices[i].mz)] = params[i];
			maxmr = std::max<int>(maxmr, math::abs(indices[i].mr));
			maxmz = std::max<int>(maxmz, math::abs(indices[i].mz));
		}
		/*if(maxmz==0) {  // dealing with the case Jz==0 -- add only two elements in m_r
			newIndices.push_back(GenFncIndex(maxmr+1, 0, 0));
			newIndices.push_back(GenFncIndex(maxmr+2, 0, 0));
			return newIndices;
		}*/

		// 2. determine the largest amplitude of coefs that are at the boundary of existing values
		double maxval = 0;
		double A1 = 0;
		int indexz = 0;
		//note terms that are not themselves off end, but have a neighbour that is
		for (int ir = 0; ir <= maxmr + 2; ir++)
			for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end() &&
					((iz <= 0 && indPairs.find(std::make_pair(ir, iz - 2)) == indPairs.end()) ||
						(iz >= 0 && indPairs.find(std::make_pair(ir, iz + 2)) == indPairs.end()) ||
						indPairs.find(std::make_pair(ir + 1, iz)) == indPairs.end())) {
					bool y = true;
					for (int k = 0;k < fracs.size();k++) {
						if (fracs[k].mz == iz) {
							y = false;
						}
					}
					if (y) {
						maxval = fmax(fabs(indPairs[std::make_pair(ir, iz)]), maxval);
						indexz = iz;
					}
				}
			}
		//printf("Index %d maxval %g\n", indexz, maxval);

		// 4. add more terms adjacent to the existing ones at the boundary, if they are large enough
		double thresh = maxval * 0.1;
		for (int ir = 0; ir <= maxmr + 3; ir++)
			for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end()
					|| (ir == 0 && iz >= 0)) continue;  // already exists or not required
				if (absvalue(indPairs, ir - 3, iz) >= thresh ||
					absvalue(indPairs, ir - 2, iz) >= thresh ||
					absvalue(indPairs, ir - 1, iz) >= thresh ||
					absvalue(indPairs, ir, iz - 2) >= thresh ||
					absvalue(indPairs, ir, iz + 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz - 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz + 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz) >= thresh)
				{   
					bool y = true;
					for (int k = 0;k < fracs.size();k++) {
						if (fracs[k].mz == iz) {
							y = false;
						}
					}if (y) {
						newIndices.push_back(GenFncIndex(ir, iz, 0));
						numadd++;
					}
					// add a term if any of its neighbours are large enough
					
				}
			}
		// 4. finish up
		//printf("%d terms added\n", numadd);
		assert(numadd > 0);
		params.resize(newIndices.size() + 2 * fracs.size());
		int last = params.size() - 1;
		for (int i = 0; i < 2 * fracs.size(); i++) {//move fracs parameters to end
			params[last - i] = params[last - i - numadd];
			params[last - i - numadd] = 0;
		}
		return newIndices;
	}
}  // namespace actions
