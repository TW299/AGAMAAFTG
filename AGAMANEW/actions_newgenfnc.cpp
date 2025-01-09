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
			AngleFinder(const GenFncIndices& _indices, const GenFncDerivs& _dSby, const Angles& _ang) :
				indices(_indices), dSby(_dSby), ang(_ang) {
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
		private:
			const GenFncIndices& indices; ///< indices of terms in generating function
			const GenFncDerivs& dSby;     ///< amplitudes of derivatives dS/dJ_{r,z,phi}
			const Angles ang;             ///< true angles
		};

	} // internal namespace
	void GenFnc::write(FILE* ofile) const {
		fprintf(ofile, "%zd\n", indices.size());
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g %g %g %g\n",
				indices[i].mr, indices[i].mz, indices[i].mphi,
				values[i], derivs[i].Jr, derivs[i].Jz, derivs[i].Jphi);
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

	Actions GenFnc::toyJ(const Actions& J, const Angles& thetaT) const {//returns toyJ(thetaT)
		Actions JT(J);
		for (unsigned int i = 0; i < indices.size(); i++) {// compute toy actions
			double val = values[i] * cos(indices[i].mr * thetaT.thetar +
				indices[i].mz * thetaT.thetaz +
				indices[i].mphi * thetaT.thetaphi);
			JT.Jr += val * indices[i].mr;
			JT.Jz += val * indices[i].mz;
			JT.Jphi += val * indices[i].mphi;
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
		}
		return theta;
	}

	Angles GenFnc::toyA(const Angles& theta) const {//real Angs -> toy Angs
		double thetaT[3], Theta[3] = { theta.thetar, theta.thetaz, theta.thetaphi };
		AngleFinder AF(indices, derivs, theta);
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
	/*
	void GenFnc::print(void) const{
		for(int i=0; i<indices.size(); i++)
			printf("%d %d %d %g\n",indices[i].mr,indices[i].mz,indices[i].mphi,values[i]);
	}*/
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

	void GenFnc::print(void) const {
		for (int i = 0; i < values.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, values[i]);
		printf("\n");
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

	GenFncIndices GenFncFit::expand(const std::vector<double>& params)
	{   /// NOTE: here we specialize for the case of axisymmetric systems!
		assert(params.size() == numParams());
		std::map< std::pair<int, int>, double > indPairs;

		// 1. determine the extent of existing grid in (mr,mz)
		int maxmr = 0, maxmz = 0;
		for (unsigned int i = 0; i < indices.size(); i++) {
			indPairs[std::make_pair(indices[i].mr, indices[i].mz)] = params[i];
			maxmr = std::max<int>(maxmr, math::abs(indices[i].mr));
			maxmz = std::max<int>(maxmz, math::abs(indices[i].mz));
		}
		GenFncIndices newIndices = indices;
		/*if(maxmz==0) {  // dealing with the case Jz==0 -- add only two elements in m_r
			newIndices.push_back(GenFncIndex(maxmr+1, 0, 0));
			newIndices.push_back(GenFncIndex(maxmr+2, 0, 0));
			return newIndices;
		}*/

		// 2. determine the largest amplitude of coefs that are at the boundary of existing values
		double maxval = 0;
		for (int ir = 0; ir <= maxmr + 2; ir++)
			for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end() &&
					((iz <= 0 && indPairs.find(std::make_pair(ir, iz - 2)) == indPairs.end()) ||
						(iz >= 0 && indPairs.find(std::make_pair(ir, iz + 2)) == indPairs.end()) ||
						indPairs.find(std::make_pair(ir + 1, iz)) == indPairs.end()))
					maxval = fmax(fabs(indPairs[std::make_pair(ir, iz)]), maxval);
			}

		// 3. add more terms adjacent to the existing ones at the boundary, if they are large enough
		double thresh = maxval * 0.1;
		int numadd = 0;
		for (int ir = 0; ir <= maxmr + 2; ir++)
			for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end() || (ir == 0 && iz >= 0))
					continue;  // already exists or not required
				if (absvalue(indPairs, ir - 2, iz) >= thresh ||
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
		assert(numadd > 0);
		return newIndices;
	}

}  // namespace actions
