#include "actions_newtorus.h"
//#define TESTIT

/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-6;

namespace actions {

	namespace {//internal

		struct PVUVSph {
			double u, v, phi, pu, pv, pphi;
		};
		struct DerivActUVSph {
			PVUVSph dbyJr, dbyJz, dbyJphi;
		};
		double** dmatrix(int n1, int n2) {
			double** f = new double* [n1];
			for (int i = 0;i < n1;i++)
				f[i] = new double[n2];
			return f;
		}
		void delmatrix(double** f, int n1) {
			for (int i = 0;i < n1;i++) delete[] f[i];
			delete[] f;
		}
		//Two functions used in containsPoint
		double angleDiff(Angles A1, Angles A2) {

			double dr = A1.thetar - A2.thetar;
			while (dr > M_PI) dr -= 2 * M_PI;
			while (dr < -M_PI) dr += 2 * M_PI;
			double dz = A1.thetaz - A2.thetaz;
			while (dz > M_PI) dz -= 2 * M_PI;
			while (dz < -M_PI) dz += 2 * M_PI;
			double dp = A1.thetaphi - A2.thetaphi;
			while (dp > M_PI) dp -= 2 * M_PI;
			while (dp < -M_PI) dp += 2 * M_PI;
			return sqrt(dr * dr + dz * dz + dp * dp);
		}

		bool is_new(Angles A1, std::vector<Angles> As) {
			bool ok = true; const double tiny = 1e-5;
			for (int i = 0; i < As.size(); i++) {
				double diff = angleDiff(A1, As[i]);
				ok = ok && diff > tiny;
			}
			return ok;
		}
		class locFinder : public math::IFunctionNdimDeriv {
		private:
			const Torus T;
			const coord::PosCyl P0;
		public:
			locFinder(const Torus _T, const coord::PosCyl& _P0) : T(_T), P0(_P0) {}
			virtual void evalDeriv(const double params[], double* values, double* derivs = NULL) const {
				Angles A1(params[0], params[1], 0);
				coord::PosCyl P1;  DerivAng<coord::Cyl> dA;
				if (derivs) {
					P1 = (T.PosDerivs(A1, dA));
					derivs[0] = dA.dbythetar.R; derivs[1] = dA.dbythetaz.R;
					derivs[2] = dA.dbythetar.z; derivs[3] = dA.dbythetaz.z;
				}
				if (values && !derivs)
					P1 = T.from_toy(A1);
				if (values) {
					values[0] = P1.R - P0.R; values[1] = P1.z - P0.z;
					double dist = sqrt(pow_2(values[0]) + pow_2(values[1]));
				}
			}
			virtual unsigned int numVars() const {
				return 2;
			}
			virtual unsigned int numValues() const {
				return 2;
			}
			void print(void) const {
				printf("LF: %f %f %f", T.Is.Js, T.Is.b, T.cs.Delta);
			}
		};

		double insertLine(int& ntop, const int tmax, double s, GenFncIndex I,
			std::vector<double>& Hmods, GenFncIndices& Hindices) {
			if (ntop == 0) {
				Hmods.push_back(s);
				Hindices.push_back(I);
				ntop++; return Hmods[0];
			}
			int l = 0;
			while (l < ntop && s <= Hmods[l]) l++;
			if (l == ntop) {//line isn't stronger than any previous line
				if (ntop >= tmax) {
					return Hmods[ntop - 1];//no room for this line
				}
				else {//add line to end of list
					Hmods.push_back(s);
					Hindices.push_back(I);
					ntop++;
					return Hmods.back();
				}
			}
			else {//we should insert current line
				if (ntop < tmax) {//move existing terms down
					Hmods.push_back(Hmods[ntop - 1]);
					Hindices.push_back(Hindices[ntop - 1]);
				}
				for (int m = ntop - 1; m > l; m--) {
					Hmods[m] = Hmods[m - 1];
					Hindices[m] = Hindices[m - 1];
				}//we've created a space at l, so fill it
				Hmods[l] = s;
				Hindices[l] = I;
				if (ntop < tmax) ntop++;//we've added rather than replaced a line
				return Hmods.back();//return weakest retained line
			}
		}
		/* compute the best focal distance at a 2d grid in L, Xi=Jz/L
		 * on input gridL, which is a grid in Jcirc, times gridXi is a uniform grid on (0,1)
		*/
		void createGridFocalDistance(const potential::BasePotential& pot,
			std::vector<double>& gridL, std::vector<double>& gridXi,
			math::Matrix<double>& grid2dD, math::Matrix<double>& grid2dR)
		{
			int sizeL = gridL.size(), sizeXi = gridXi.size();
			math::Matrix<double> grid2dL(sizeL, sizeXi);
			for (int iL = 1; iL < sizeL - 1; iL++) {//omit bdy values
				double Rc, Vc, Jz, Jc = gridL[iL];
				double E = E_circ(pot, Jc, &Rc, &Vc);
				std::vector<double> L_vals(sizeXi);
				std::vector<double> Xi_vals(sizeXi);
				std::vector<double> D_vals(sizeXi);
				std::vector<double> R_vals(sizeXi);
				for (int iXi = 1; iXi < sizeXi - 1; iXi++) {
					double Jphi = Jc * (1 - gridXi[iXi]);
					double Rsh, FD;
					FD = estimateFocalDistanceShellOrbit(pot, E, Jphi, &Rsh, &Jz);
					double L = fabs(Jphi) + Jz;
					L_vals[iXi] = L; Xi_vals[iXi] = Jz / L;
					D_vals[iXi] = FD; R_vals[iXi] = Rsh;
					if (iXi > 0 && Xi_vals[iXi] < Xi_vals[iXi - 1])
						printf("createGridFocalDistances: non-monotonic Xi:\n",
							"%d %f %f\n", iXi, Xi_vals[iXi - 1], Xi_vals[iXi]);
				}
				// bdy values
				L_vals[0] = gridL[iL]; L_vals[sizeXi - 1] = L_vals[sizeXi - 2];
				Xi_vals[0] = 0; Xi_vals[sizeXi - 1] = 1;
				D_vals[0] = D_vals[1]; D_vals[sizeXi - 1] = D_vals[sizeXi - 2];
				R_vals[0] = R_vals[1]; R_vals[sizeXi - 1] = R_vals[sizeXi - 2];
				math::LinearInterpolator interpL(Xi_vals, L_vals);
				math::LinearInterpolator interpD(Xi_vals, D_vals);
				math::LinearInterpolator interpR(Xi_vals, R_vals);
				for (int iXi = 0; iXi < sizeXi; iXi++) {//interpolae L, D onto regular grid in Xi
					interpL.evalDeriv(gridXi[iXi], &grid2dL(iL, iXi));
					interpD.evalDeriv(gridXi[iXi], &grid2dD(iL, iXi));
					interpR.evalDeriv(gridXi[iXi], &grid2dR(iL, iXi));
				}
			}
			//We now need to fill in rows iL=0, iL=sizeL-1
			for (int iXi = 0; iXi < sizeXi; iXi++) {
				grid2dL(0, iXi) = 0; grid2dL(sizeL - 1, iXi) = 1.01 * grid2dL(sizeL - 2, iXi);
				grid2dD(0, iXi) = 0; grid2dD(sizeL - 1, iXi) = grid2dD(sizeL - 2, iXi);
				grid2dR(0, iXi) = 0; grid2dR(sizeL - 1, iXi) = grid2dR(sizeL - 2, iXi);
			}
			//Now grid2dD contains D on regular grid in Xi but irregular
			//values of L that are stored in grid2dL
			for (int iXi = 0; iXi < sizeXi; iXi++) {
				std::vector<double> L_vals(sizeL);
				std::vector<double> D_vals(sizeL);
				std::vector<double> R_vals(sizeL);
				for (int iL = 0; iL < sizeL; iL++) {
					L_vals[iL] = grid2dL(iL, iXi);
					D_vals[iL] = grid2dD(iL, iXi);
					R_vals[iL] = grid2dR(iL, iXi);
					if (iL > 0 && L_vals[iL] <= L_vals[iL - 1])
						printf("createGridFocalDistance: non-monotonic L_vals %g %g\n",
							L_vals[iL - 1], L_vals[iL]);
				}
				math::LinearInterpolator DL(L_vals, D_vals);
				math::LinearInterpolator RL(L_vals, R_vals);
				for (int iL = 0; iL < sizeL; iL++) {
					DL.evalDeriv(gridL[iL], &grid2dD(iL, iXi));
					RL.evalDeriv(gridL[iL], &grid2dR(iL, iXi));
					//if(std::isnan(grid2dD(iL,iXi))) printf("Dgrid nan %d %d\n",iL,iXi);
					//if(std::isnan(grid2dR(iL,iXi))) printf("Rgrid nan %d %d\n",iL,iXi);
				}
			}
		}

		double H_dHdX(const potential::BasePotential& pot, const coord::PosMomCyl Rzphi,
			coord::PosMomCyl& dHdX) {
			double Phi; coord::GradCyl grad;
			pot.eval(Rzphi, &Phi, &grad);
			dHdX.R = grad.dR - pow_2(Rzphi.pphi / Rzphi.R) / Rzphi.R;
			dHdX.z = grad.dz; dHdX.phi = grad.dphi;
			dHdX.pR = Rzphi.pR; dHdX.pz = Rzphi.pz; dHdX.pphi = Rzphi.pphi / pow_2(Rzphi.R);
			return .5 * (pow_2(Rzphi.pR) + pow_2(Rzphi.pz) + pow_2(Rzphi.pphi / Rzphi.R)) + Phi;
		}

		coord::PosMomCyl PosMomDerivs(const coord::PosMomSph& rtheta, const coord::UVSph& cs,
			const actions::DerivAct<coord::Sph>* dJ,//input from Is
			actions::DerivAct<coord::Cyl>& dJC) {
			double snv, csv;
			math::sincos(rtheta.theta, snv, csv);
			coord::PosUVSph uv(asinh(rtheta.r / cs.Delta), rtheta.theta, rtheta.phi, cs);
			double chu = cosh(uv.u), shu = sinh(uv.u);
			double dudr = 1 / (cs.Delta * chu);
			coord::MomUVSph pp(rtheta.pr / dudr, rtheta.ptheta, rtheta.pphi);
			coord::PosMomUVSph uvp(uv, pp);
			coord::PosMomCyl Rzphi = coord::toPosMomCyl(uvp);
			if (!dJ) return Rzphi;
			double dpudr = cs.Delta * shu * dudr * rtheta.pr;
			double dpvdtheta = 0;
			double dpudpr = 1 / dudr;
			double dpvdptheta = 1;
			DerivActUVSph dJUV;
			dJUV.dbyJr.u = dudr * dJ->dbyJr.r;
			dJUV.dbyJz.u = dudr * dJ->dbyJz.r;
			dJUV.dbyJphi.u = dudr * dJ->dbyJphi.r;
			dJUV.dbyJr.v = dJ->dbyJr.theta;
			dJUV.dbyJz.v = dJ->dbyJz.theta;
			dJUV.dbyJphi.v = dJ->dbyJphi.theta;
			dJUV.dbyJr.phi = dJ->dbyJr.phi;//all 3 non-zero
			dJUV.dbyJz.phi = dJ->dbyJz.phi;
			dJUV.dbyJphi.phi = dJ->dbyJphi.phi;

			dJUV.dbyJr.pu = dpudr * dJ->dbyJr.r + dpudpr * dJ->dbyJr.pr;
			dJUV.dbyJz.pu = dpudr * dJ->dbyJz.r + dpudpr * dJ->dbyJz.pr;
			dJUV.dbyJphi.pu = dpudr * dJ->dbyJphi.r + dpudpr * dJ->dbyJphi.pr;
			dJUV.dbyJr.pv = dJ->dbyJr.ptheta;
			dJUV.dbyJz.pv = dJ->dbyJz.ptheta;
			dJUV.dbyJphi.pv = dJ->dbyJphi.ptheta;
			dJUV.dbyJr.pphi = dJ->dbyJr.pphi;//both these vanish!
			dJUV.dbyJz.pphi = dJ->dbyJz.pphi;
			dJUV.dbyJphi.pphi = dJ->dbyJphi.pphi;//unity

			double dRdu = cs.Delta * chu * snv, dRdv = cs.Delta * shu * csv;
			double dzdu = cs.Delta * shu * csv, dzdv = -cs.Delta * chu * snv;
			double det = cs.Delta * (pow_2(shu) + pow_2(snv));
			double ddetdu = 2 * shu * chu * cs.Delta, ddetdv = 2 * snv * csv * cs.Delta;
			double dpRdu = (shu * snv * uvp.pu + chu * csv * uvp.pv - Rzphi.pR * ddetdu) / det;
			double dpRdv = (chu * csv * uvp.pu - shu * snv * uvp.pv - Rzphi.pR * ddetdv) / det;
			double dpzdu = (chu * csv * uvp.pu - shu * snv * uvp.pv - Rzphi.pz * ddetdu) / det;
			double dpzdv = (-shu * snv * uvp.pu - chu * csv * uvp.pv - Rzphi.pz * ddetdv) / det;
			double dpRdpu = chu * snv / det, dpRdpv = shu * csv / det;
			double dpzdpu = shu * csv / det, dpzdpv = -chu * snv / det;

			dJC.dbyJr.R = dRdu * dJUV.dbyJr.u + dRdv * dJUV.dbyJr.v;
			dJC.dbyJz.R = dRdu * dJUV.dbyJz.u + dRdv * dJUV.dbyJz.v;
			dJC.dbyJphi.R = dRdu * dJUV.dbyJphi.u + dRdv * dJUV.dbyJphi.v;
			dJC.dbyJr.z = dzdu * dJUV.dbyJr.u + dzdv * dJUV.dbyJr.v;
			dJC.dbyJz.z = dzdu * dJUV.dbyJz.u + dzdv * dJUV.dbyJz.v;
			dJC.dbyJphi.z = dzdu * dJUV.dbyJphi.u + dzdv * dJUV.dbyJphi.v;
			dJC.dbyJr.phi = dJUV.dbyJr.phi;  //dphidu*dJUV.dbyJr.u  + dphidv*dJUV.dbyJr.v;
			dJC.dbyJz.phi = dJUV.dbyJz.phi;  //dphidu*dJUV.dbyJz.u  + dphidv*dJUV.dbyJz.v;
			dJC.dbyJphi.phi = dJUV.dbyJphi.phi;//dphidu*dJUV.dbyJphi.u+ dphidv*dJUV.dbyJphi.v;

			dJC.dbyJr.pR = dpRdu * dJUV.dbyJr.u + dpRdv * dJUV.dbyJr.v + dpRdpu * dJUV.dbyJr.pu + dpRdpv * dJUV.dbyJr.pv;
			dJC.dbyJz.pR = dpRdu * dJUV.dbyJz.u + dpRdv * dJUV.dbyJz.v + dpRdpu * dJUV.dbyJz.pu + dpRdpv * dJUV.dbyJz.pv;
			dJC.dbyJphi.pR = dpRdu * dJUV.dbyJphi.u + dpRdv * dJUV.dbyJphi.v + dpRdpu * dJUV.dbyJphi.pu + dpRdpv * dJUV.dbyJphi.pv;
			dJC.dbyJr.pz = dpzdu * dJUV.dbyJr.u + dpzdv * dJUV.dbyJr.v + dpzdpu * dJUV.dbyJr.pu + dpzdpv * dJUV.dbyJr.pv;
			dJC.dbyJz.pz = dpzdu * dJUV.dbyJz.u + dpzdv * dJUV.dbyJz.v + dpzdpu * dJUV.dbyJz.pu + dpzdpv * dJUV.dbyJz.pv;
			dJC.dbyJphi.pz = dpzdu * dJUV.dbyJphi.u + dpzdv * dJUV.dbyJphi.v + dpzdpu * dJUV.dbyJphi.pu + dpzdpv * dJUV.dbyJphi.pv;
			dJC.dbyJr.pphi = 0;
			dJC.dbyJz.pphi = 0;
			dJC.dbyJphi.pphi = 1;
			if (std::isnan(dJC.dbyJr.pR)) {
				printf("PosMomD\n%f %f %f %f %f %f %f %f\n", dpRdu, dJUV.dbyJr.u, dpRdv, dJUV.dbyJr.v,
					dpRdpu, dJUV.dbyJr.pu, dpRdpv, dJUV.dbyJr.pv);
				exit(0);
			}
			return Rzphi;
		}

		class Iso {
		private:
			double L, Jr;
		public:
			Iso(double _L, double _Jr) : L(_L), Jr(_Jr) {}
			double b2cE(double Js) const {//general E
				return .5 * pow_2(Js * Js) / pow_2(Jr + .5 * (L + sqrt(L * L + 4 * Js * Js)));
			}
			double b2cEc(double Js) const {//circular E
				return .5 * pow_2(Js * Js) / pow_2(.5 * (L + sqrt(L * L + 4 * Js * Js)));
			}
			double g(double Js) const {//Rsh=b*g(Js)
				double g2 = (pow((L + sqrt(L * L + 4 * Js * Js)) / (2 * Js), 4) - 1);
				if (g2 < 0) printf("g2<0: %f\n", g2);
				return sqrt(g2);
			}
			double cob(double Js) const {//c/b
				return .5 * pow_2(Js) / b2cE(Js) - 1;;
			}
			double ecc(double Js) const {
				double boc = 1 / cob(Js);
				double e2 = 1 - pow_2(L / Js) * boc * (1 + boc);
				return e2 < 1 ? sqrt(e2) : 0;
			}
			double f(double b, double Js, double& e) const {//ratio of forces aopo/peri
				double cb = cob(Js), boc = 1 / cb;
				e = ecc(Js);
				double up = 1 + e, um = 1 - e;
				double rp = b * cb * sqrt(up * (up + 2 * boc)), ap = sqrt(b * b + rp * rp);
				double rm = b * cb * sqrt(um * (um + 2 * boc)), am = sqrt(b * b + rm * rm);
				return pow_2((b + am) / (b + ap)) * am / ap * rp / rm;
			}
		};
		/*
		 Class to find Js (and thus b) of isochrone for which the isochrone's
		 force ratio matches that in real Phi
		*/
		class JsFinder : public math::IFunction {
			const potential::BasePotential& pot;
			const Iso Is;
			const double Rsh;
		public:
			JsFinder(const potential::BasePotential& _pot, const Iso& _Is, const double _Rsh) :
				pot(_pot), Is(_Is), Rsh(_Rsh) {
			}
			virtual void evalDeriv(double Js, double* val, double* deriv = 0, double* deriv2 = 0) const {
				double b = Rsh / Is.g(Js);
				double e, f_apo_peri = Is.f(b, Js, e);
				double c = Is.cob(Js) * b;
				double F[2];
				for (int k = -1; k < 2; k += 2) {
					double u = 1 + k * e;
					double Phi, r = c * sqrt(u * (u + 2 * b / c));
					coord::PosCyl Rz(r, 0, 0); coord::GradCyl grad;
					pot.eval(Rz, &Phi, &grad);
					F[(k + 1) / 2] = grad.dR;
				}
				*val = f_apo_peri - F[1] / F[0];
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}

		};

		/// create the array of indices of the generating function with all terms up to the given maximum order
		static GenFncIndices makeGridIndices(int irmax, int izmax)
		{   /// NOTE: here we specialize for the case of axisymmetric systems!
			GenFncIndices indices;
			for (int ir = 0; ir <= irmax; ir++)
				for (int iz = -izmax; iz <= (ir == 0 ? -2 : izmax); iz += 2) {
					indices.push_back(GenFncIndex(ir, iz, 0));
				}
					
					
			return indices;
		}

		/// return the absolute value of an element in a map, or zero if it doesn't exist
		static inline double absvalue(const std::map< std::pair<int, int>, double >& indPairs, int ir, int iz)
		{
			if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
				return fabs(indPairs.find(std::make_pair(ir, iz))->second);
			else
				return 0;
		}
		/** Compute the derivative of Hamiltonian by toy actions:
			dH/dJ = dH/d{x,v} d{x,v}/dJ, where the lhs is a covector of length 3,
			the first term on rhs is a covector of length 6 (the gradient dPhi/dx and the velocity),
			and the second term is a 6x3 matrix of partial derivs provided by the toy map.
		*/
		static inline Actions dHbydJ(const coord::PosMomCyl& dHby,
			const DerivAct<coord::Cyl>& dXdJ)
		{
			return Actions(
				dXdJ.dbyJr.R * dHby.R + dXdJ.dbyJr.z * dHby.z + dXdJ.dbyJr.phi * dHby.phi +
				dXdJ.dbyJr.pR * dHby.pR + dXdJ.dbyJr.pz * dHby.pz + dXdJ.dbyJr.pphi * dHby.pphi,
				dXdJ.dbyJz.R * dHby.R + dXdJ.dbyJz.z * dHby.z + dXdJ.dbyJz.phi * dHby.phi +
				dXdJ.dbyJz.pR * dHby.pR + dXdJ.dbyJz.pz * dHby.pz + dXdJ.dbyJz.pphi * dHby.pphi,
				dXdJ.dbyJphi.R * dHby.R + dXdJ.dbyJphi.z * dHby.z + dXdJ.dbyJphi.phi * dHby.phi +
				dXdJ.dbyJphi.pR * dHby.pR + dXdJ.dbyJphi.pz * dHby.pz + dXdJ.dbyJphi.pphi * dHby.pphi);
		}
		//class to help with SoS computation
		class CrossingFinder : public math::IFunction {
		private:
			const Torus* T;
			const double thetaT_r;
		public:
			CrossingFinder(const Torus* _T, const double& _thetaT) : T(_T), thetaT_r(_thetaT) {}
			virtual void evalDeriv(const double thetaT_z, double* value, double* deriv = NULL, double* deriv2 = NULL) const {
				Angles thetaT(thetaT_r, thetaT_z, 0);
				coord::PosMomCyl Rz(T->from_toy(thetaT));
				*value = Rz.z;
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}
		};
		/*
		 * Function to comput matrix dJ/dp = dx/dthetaT. It
		 * takes in dx/dthataT and dtheta/dthetaT, inverts the
		 * latter and multiplies on the former
		 */
		void assemble(DerivAng<coord::Cyl>& dA, math::Matrix<double>& M) {
			std::vector<double> dxdthetaT = { dA.dbythetar.R,
			dA.dbythetaz.R,dA.dbythetaphi.R,dA.dbythetar.z,
			dA.dbythetaz.z,dA.dbythetaphi.z,dA.dbythetar.phi,
			dA.dbythetaz.phi,dA.dbythetaphi.phi };
			math::Matrix<double> dxdthetaT2(3, 3);
			for (int k = 0;k < 3;k++) {
				for (int n = 0;n < 3;n++) {
					dxdthetaT2(k, n) = dxdthetaT[3 * k + n];
				}
			}
			math::LUDecomp L2(M);
			math::Matrix<double> M2 = (L2.inverse(3));
			math::Matrix<double> Mat3(3, 3);
			math::CBLAS_TRANSPOSE Transp = math::CblasNoTrans;
			math::blas_dgemm(math::CblasNoTrans, math::CblasNoTrans, 1.0, dxdthetaT2, M2, 0.0, Mat3);
			//DerivAng<coord::Cyl> deriv1;

			dA.dbythetar.R = Mat3.at(0, 0);dA.dbythetaz.R = Mat3.at(0, 1);dA.dbythetaphi.R = Mat3.at(0, 2);
			dA.dbythetar.z = Mat3.at(1, 0);dA.dbythetaz.z = Mat3.at(1, 1);dA.dbythetaphi.z = Mat3.at(1, 2);
			dA.dbythetar.phi = Mat3.at(2, 0);dA.dbythetaz.phi = Mat3.at(2, 1);dA.dbythetaphi.phi = Mat3.at(2, 2);
			//printf("dA %f %f\n",dA.dbythetar.R,dA.dbythetar.z);
		}
		/* A helper class used by TorusGenerator. It evaluates H and its
		 * derivatives on a toy-angle grid over a proposed torus using a
		 * proposed GF. Once the GF has been
		 * optimised, it runs over the grid a final time to determine d S_k/d J
		 * so we can recover the true angles
		 */
		class torusFitter : public math::IFunctionNdimDeriv {
		private:
			const Actions& J;
			const potential::BasePotential& pot;
			const Isochrone& Is;
			const coord::UVSph& cs;
			GenFncFitSeries& GFFS;
			double NANfrac;
		public:
			torusFitter(const Actions& _J,
				const potential::BasePotential& _pot,
				const Isochrone& _Is, const coord::UVSph& _cs,
				GenFncFitSeries& _GFFS) :
				J(_J), pot(_pot), Is(_Is), cs(_cs), GFFS(_GFFS), NANfrac(0) {
			}
			virtual unsigned int numVars() const {
				return GFFS.numParams();
			}
			virtual unsigned int numValues() const {
				return GFFS.numPoints();
			}
			double computeHamiltonianAtPoint(const double params[],
				const unsigned int indPoint, Actions* dHdJ = NULL,
				double* derivGenFnc = NULL) const
			{
				// Generating function computes the toy actions from the real actions
				// at the given point in the grid of toy angles grid
				ActionAngles toyAA = GFFS.toyActionAngles(indPoint, params);

				// do not allow to stray into forbidden region of negative actions
					//assert((toyAA.Jr >= 0) && (toyAA.Jz >= 0),"negative action in computeHamiltonianAtPoint");
				if (toyAA.Jr < 0 || toyAA.Jz < 0) {
					//printf("J<0: %f %f ", toyAA.Jr, toyAA.Jz);
					return NAN;
				}
				// aa2pq computes the position and velocity from the toy actions and angles,
				// and optionally their derivatives w.r.t. toy actions and toy map parameters,
				DerivAct<coord::Sph> derivAct;//derivs of coords wrt Js
				coord::PosMomSph rtheta = Is.aa2pq(toyAA, NULL, derivGenFnc != NULL ? &derivAct : NULL);
				actions::DerivAct<coord::Cyl> dXdJ;
				coord::PosMomCyl Rzphi = PosMomDerivs(rtheta, cs,
					derivGenFnc != NULL ? &derivAct : NULL, dXdJ);
				// obtain the value of the real Hamiltonian at the given point and its
				// derivatives w.r.t. coordinates/momenta
				coord::PosMomCyl dHdX;
				double H = H_dHdX(pot, Rzphi, dHdX);

				// derivatives of Hamiltonian w.r.t. parameters of gen.fnc.
				if (derivGenFnc) {
					Actions dHby = dHbydJ(dHdX, dXdJ);// derivative of Hamiltonian by toy actions
					if (indPoint == -30 || indPoint == -10) {
						printf("(H %f dHby %f %f\n", H, dHby.Jr, dHby.Jz);
						printf("dXdJ %f %f %f %f %f %f %f %f\n", dXdJ.dbyJr.R, dXdJ.dbyJr.z, dXdJ.dbyJz.R, dXdJ.dbyJz.z, dXdJ.dbyJr.pR, dXdJ.dbyJr.pz, dXdJ.dbyJz.pR, dXdJ.dbyJz.pz);
						printf("[dHdX %f %f %f %f %f %f\n] ", dHdX.R, dHdX.z, dHdX.phi, dHdX.pR, dHdX.pz, dHdX.pphi);
					}
					if (dHdJ) *dHdJ = dHby;
					for (unsigned int p = 0; p < GFFS.numParams(); p++) {
						// derivs of toy actions by gen.fnc.params
						Actions dbyS = GFFS.deriv(indPoint, p, &params[0]);
						// derivs of Hamiltonian by gen.fnc.params
						double  dHdS = dHby.Jr * dbyS.Jr + dHby.Jz * dbyS.Jz + dHby.Jphi * dbyS.Jphi;
						derivGenFnc[p] = dHdS;
						if (std::isnan(dHdS)) {
							printf("AtPoint  %f %f %f %f %f %f\n", dHby.Jr, dHby.Jz, dHby.Jphi, dbyS.Jr, dbyS.Jz, dbyS.Jphi);
							printf("%f %f %f %f %f %f\n", dXdJ.dbyJr.R, dXdJ.dbyJr.z, dXdJ.dbyJr.phi,
								dXdJ.dbyJr.pR, dXdJ.dbyJr.pz, dXdJ.dbyJr.pphi);
							printf("%f %f %f %f %f %f\n", dXdJ.dbyJz.R, dXdJ.dbyJz.z, dXdJ.dbyJz.phi,
								dXdJ.dbyJz.pR, dXdJ.dbyJz.pz, dXdJ.dbyJz.pphi);
							printf("%f %f %f %f %f %f\n", dXdJ.dbyJphi.R, dXdJ.dbyJphi.z, dXdJ.dbyJphi.phi,
								dXdJ.dbyJphi.pR, dXdJ.dbyJphi.pz, dXdJ.dbyJphi.pphi);
							exit(0);
						}
					}
				}
				return H;
			}

			double computeHamiltonianDisp(const std::vector<double>& params, double& Hbar)
			{
				const unsigned int numParams = GFFS.numParams(), max_threads = 1;
				double Hsm[max_threads] = { 0 }, Hsq[max_threads] = { 0 };
				int N[max_threads] = { 0 };
				int Nnan[max_threads] = { 0 };
//#pragma omp parallel for schedule(dynamic)
				for (int indPoint = 0; indPoint < GFFS.numPoints(); indPoint++) {
					//int nth = omp_get_thread_num();
					double H = computeHamiltonianAtPoint(&params[0], indPoint);
					if (std::isnan(H)) Nnan[0]++;
					else {
						//Hsm[nth] += H; Hsq[nth] += H * H; N[nth]++;
						Hsm[0] += H; Hsq[0] += H * H; N[0]++;
					}
				}
				/*for (int i = 1; i < max_threads; i++) {//concatenate sums
					Nnan[0] += Nnan[i]; Hsm[0] += Hsm[i]; Hsq[0] += Hsq[i]; N[0] += N[i];
				}*/
				Hbar = Hsm[0] / N[0];
				Hsq[0] = Hsq[0] / N[0] - Hbar * Hbar;
				NANfrac = (double)Nnan[0] / (double)GFFS.numPoints();
				return (Hsq[0]);
			}
			double giveNANfrac() const {
				return NANfrac;
			}
			void evalDeriv(const double params[],
				double* deltaHvalues, double* dHdParams) const
			{
				const unsigned int numPoints = GFFS.numPoints();
				const unsigned int numParams = GFFS.numParams();
				const unsigned int max_threads = 1;

				// we need to store the values of Hamiltonian at grid points even if this is not requested,
				// because they are used to correct the entries of the Jacobian matrix
				// to account for the fact that the mean <H> also depends on the parameters
				std::vector<double> Hvalues(numPoints);
				double Havg[max_threads] = { 0 };  // accumulator for the average Hamiltonian
				int nNAN[max_threads] = { 0 };
				// loop over grid of toy angles
//#pragma omp parallel for schedule(dynamic)
				for (int indPoint = 0; indPoint < numPoints; indPoint++) {
					//int nth = omp_get_thread_num();
					double H = computeHamiltonianAtPoint(params, indPoint, NULL,
						dHdParams ? dHdParams + indPoint * numParams : NULL);
					if (std::isnan(H)) {
						/*nNAN[nth]++;
						Havg[nth] += 1e6;
						Hvalues[indPoint] = 1e6; break;
						//*/
						nNAN[0]++;
						Havg[0] += 1e6;
						Hvalues[indPoint] = 1e6; break;
					}
					else {
						// accumulate the average value and store the output
						//Havg[nth] += H;
						Havg[0] += H;
						Hvalues[indPoint] = H;
					}
				}
				/*for (int i = 1; i < max_threads; i++) {
					nNAN[0] += nNAN[i]; Havg[0] += Havg[i];
				}*/

				// convert from  H_k  to  deltaH_k = H_k - <H>
				nNAN[0] /= numPoints; Havg[0] /= numPoints;//we can't put nNAN in NANfrac!
				if (deltaHvalues) {
					double disp = 0;
					for (unsigned int indPoint = 0; indPoint < numPoints; indPoint++) {
						deltaHvalues[indPoint] = Hvalues[indPoint] - Havg[0];
						disp += pow_2(deltaHvalues[indPoint]);
					}
				}
				// convert derivatives:  d(deltaH_k) / dP_p = dH_k / dP_p - d<H> / dP_p
				if (dHdParams) {
					std::vector<double> dHavgdP(numPoints);
					for (unsigned int p = 0; p < numParams; p++) dHavgdP[p] = 0;
					for (unsigned int pp = 0; pp < numPoints * numParams; pp++) {
						dHavgdP[pp % numParams] += dHdParams[pp] / numPoints;
						if (std::isnan(dHavgdP[pp % numParams])) {
							printf("nan@ %d %d %g\n", pp % numParams, pp / numParams, dHdParams[pp]);
							exit(0);
						}
					}
					for (unsigned int pp = 0; pp < numPoints * numParams; pp++) {
						unsigned int indPoint = pp / numParams;
						unsigned int indParam = pp % numParams;
						dHdParams[pp] = dHdParams[pp] - dHavgdP[indParam];
					}
				}
			}
			/** Compute the frequencies and the derivatives of generating function by real actions,
				used in angle mapping.
				The three arrays of derivatives dS_i/dJ_{r,z,phi}, i=0..numParamsGenFnc-1,
				together with three frequencies Omega_{r,z,phi}, are the solutions of
				an overdetermined system of linear equations:
				\f$  M_{k,i} X_{i} = RHS_{k}, k=0..numPoints-1  \f$,
				where numPoints is the number of individual triplets of toy angles,
				\f$  X_i  \f$ is the solution vector {Omega, dS_i/dJ} for each direction (r,z,theta),
				\f$  RHS_k = dH/dJ(\theta_k)  \f$, again for three directions independently, and
				\f$  M_{k,i}  \f$ is the matrix of coefficients shared between all three equation systems:
				\f$  M_{k,0} = 1, M_{k,i+1} = -dH/dS_i(\theta_k)  \f$.
				The matrix M and three RHS vectors are filled using the same approach as during
				the Levenberg-Marquardt minimization, from the provided parameters of toy map and
				generating function; then the three linear systems are solved using
				the singular-value decomposition of the shared coefficient matrix,
				and the output frequencies and gen.fnc.derivatives are returned in corresponding arguments.
				The return value of this function is the same as `computeHamiltonianDisp()`.
			*/
			double fitAngleMap(const double params[],
				double& Hbar, Frequencies& freqs, GenFncDerivs& derivs) const {
				unsigned int numPoints = GFFS.numPoints();
				unsigned int numParams = GFFS.numParams();
				// the matrix of coefficients shared between three linear systems
				math::Matrix<double> coefsdHdS(numPoints, numParams + 1);
				// tmp storage for dH/dS
				std::vector<double> derivGenFnc(numParams + 1);
				// derivs of Hamiltonian by toy actions (RHS vectors)
				std::vector<double> dHdJr(numPoints), dHdJz(numPoints), dHdJphi(numPoints);
				// accumulator for computing dispersion in H
				math::Averager Havg;

				// loop over grid of toy angles
				for (unsigned int indPoint = 0; indPoint < numPoints; indPoint++) {
					Actions dHby;  // derivative of Hamiltonian by toy actions
					double H = computeHamiltonianAtPoint(&params[0], indPoint,
						&dHby, &derivGenFnc.front());
					if (!std::isnan(H)) Havg.add(H);
					else printf("nan @ %d ", indPoint);
					// fill the elements of each of three rhs vectors
					dHdJr[indPoint] = dHby.Jr;
					dHdJz[indPoint] = dHby.Jz;
					dHdJphi[indPoint] = dHby.Jphi;
					// fill the matrix row
					coefsdHdS(indPoint, 0) = 1;  // matrix coef for omega
					for (unsigned int p = 0; p < numParams; p++)
						coefsdHdS(indPoint, p + 1) = -derivGenFnc[p];  // matrix coef for dS_p/dJ
				}
				Hbar = Havg.mean();
				// solve the overdetermined linear system in the least-square sense:
				// step 1: prepare the SVD of coefs matrix
				math::SVDecomp SVD(coefsdHdS);

				// step 2: solve three linear systems with the same matrix but different rhs
				std::vector<double> dSdJr(SVD.solve(dHdJr)), dSdJz(SVD.solve(dHdJz)),
					dSdJphi(SVD.solve(dHdJphi));

				// store output
				freqs.Omegar = dSdJr[0];
				freqs.Omegaz = dSdJz[0];
				freqs.Omegaphi = dSdJphi[0];
				derivs.resize(numParams);
				for (unsigned int p = 0; p < numParams; p++) {
					derivs[p].Jr = dSdJr[p + 1];
					derivs[p].Jz = dSdJz[p + 1];
					derivs[p].Jphi = dSdJphi[p + 1];
				}
				return sqrt(Havg.disp());
			}

		};
	}//internal

	coord::PosMomCar xyPointTrans::rp2xp(const coord::PosMomSph rp) const {
		double spsi, cpsi, dpsidphi; math::sincos(PC.Psi(rp.phi, dpsidphi), spsi, cpsi);
		double x = rp.r * cpsi, y = sqrt(pow_2(rp.r) + Delta2) * spsi;
		double r2 = pow_2(rp.r);
		double A = (r2 + Delta2 * pow_2(cpsi)) * dpsidphi, rt = sqrt(r2 + Delta2);
		double px = ((r2 + Delta2) * cpsi * dpsidphi * rp.pr - rp.r * spsi * rp.pphi) / A;
		double py = rt * (rp.r * spsi * dpsidphi * rp.pr + cpsi * rp.pphi) / A;
		return coord::PosMomCar(x, y, 0, px, py, 0);
	}
	coord::PosMomSph xyPointTrans::xp2rp(const coord::PosMomCar xp) const {
		double x2 = pow_2(xp.x), y2 = pow_2(xp.y), B = x2 + y2 - Delta2;
		double r2 = .5 * (B + sqrt(B * B + 4 * x2 * Delta2)), r = sqrt(r2);
		double cpsi = xp.x / r, rt = sqrt(r2 + Delta2);
		double spsi = xp.y > 0 ? sqrt(1 - cpsi * cpsi) : -sqrt(1 - cpsi * cpsi);
		double pr = cpsi * xp.px + r * spsi / rt * xp.py;
		double psi = atan2(spsi, cpsi), dpsidphi, phi = PC.Phi(psi, dpsidphi);
		double pphi = (-r * spsi * xp.px + rt * cpsi * xp.py) * dpsidphi;
		return coord::PosMomSph(r, .5 * M_PI, phi, pr, 0, pphi);
	}
	coord::PosMomSph PointTrans::Cyl2Sph(const coord::PosMomCyl Rz) const {
		double R2 = pow_2(Rz.R), z2 = pow_2(Rz.z);
		double B = R2 + z2 - cs.Delta2;
		double r2 = .5 * (B + sqrt(B * B + 4 * R2 * cs.Delta2)), r = sqrt(r2);
		double rt = sqrt(r2 + cs.Delta2);
		double theta = acos(Rz.z / rt);
		double snt, cst; math::sincos(theta, snt, cst);
		double pr = snt * Rz.pR + r / rt * cst * Rz.pz;
		double ptheta = r * cst * Rz.pR - rt * snt * Rz.pz;
		return coord::PosMomSph(r, theta, Rz.phi, pr, ptheta, Rz.pphi);
	}

	coord::PosMomCyl PointTrans::Sph2Cyl(const coord::PosMomSph rp) const {
		double sq = pow_2(rp.r) + cs.Delta2, rt = sqrt(sq);
		double snt, cst; math::sincos(rp.theta, snt, cst);
		double R = rp.r * snt, z = rt * cst;
		double pR = (sq * snt * rp.pr + rp.r * cst * rp.ptheta) / (sq * pow_2(snt) + pow_2(rp.r * cst));
		double pz = (rp.r * cst * rp.pr - snt * rp.ptheta) / (pow_2(rp.r * cst) / rt + rt * pow_2(snt));
		return coord::PosMomCyl(R, z, rp.phi, pR, pz, rp.pphi);
	}

	//interpolate between 2 tori
	EXP Torus interpT(const double x, Torus T0, Torus T1) {
		if (x == 1) return T0;
		if (x == 0) return T1;
		Torus T2(T0 * x);
		T2 += T1 * (1 - x);
		return T2;
	}
	EXP Torus InterpTorus(std::vector<Torus>& Tgrid, std::vector<double>& xs, double x) {
		int top, bot;
		if (xs[0] > xs[xs.size() - 1]) {
			top = 0; bot = xs.size() - 1;
		}
		else {
			top = xs.size() - 1; bot = 0;
		}//now top should point to largest x
		double f;
		if ((xs[top] - x) * (x - xs[bot]) < 0.) {//x lies out of grid
			if (x < xs[bot]) {
				top = bot + 1; f = 1;//T = Tgrid[bot];
			}
			else {
				bot = top - 1; f = 0;// T = Tgrid[top];
			}
		}
		else {
			while (abs(top - bot) > 1) {
				int n = (top + bot) / 2;
				if ((xs[top] - x) * (x - xs[n]) >= 0) bot = n;
				else top = n;
			}
			f = (xs[top] - x) / (xs[top] - xs[bot]);//distance from top
		}
		return interpT(f, Tgrid[bot], Tgrid[top]);
	}
	EXP eTorus interpT(const double x, const eTorus T0, const eTorus T1) {
		if (x == 1) return T0;
		if (x == 0) return T1;
		eTorus T2(T0 * x);
		T2 += T1 * (1 - x);
		return T2;
	}
	EXP eTorus InterpTorus(std::vector<eTorus>& Tgrid, std::vector<double>& xs, double x) {
		int top, bot;
		if (xs[0] > xs[xs.size() - 1]) {
			top = 0; bot = xs.size() - 1;
		}
		else {
			top = xs.size() - 1; bot = 0;
		}//now top should point to largest x
		double f;
		if ((xs[top] - x) * (x - xs[bot]) < 0.) {//x lies out of grid
			if (x < xs[bot]) {
				top = bot + 1; f = 1;//T = Tgrid[bot];
			}
			else {
				bot = top - 1; f = 0;// T = Tgrid[top];
			}
		}
		else {
			while (abs(top - bot) > 1) {
				int n = (top + bot) / 2;
				if ((xs[top] - x) * (x - xs[n]) >= 0) bot = n;
				else top = n;
			}
			f = (xs[top] - x) / (xs[top] - xs[bot]);//distance from top
		}
		return interpT(f, Tgrid[bot], Tgrid[top]);
	}
	PerturbingHamiltonian& PerturbingHamiltonian::operator *= (const double a) {
		for (int i = 0; i < values.size(); i++)
			values[i] *= a;
		return *this;
	}
	PerturbingHamiltonian PerturbingHamiltonian::operator * (const double a) const {
		PerturbingHamiltonian H2(indices, values);
		H2 *= a;
		return H2;
	}
	PerturbingHamiltonian& PerturbingHamiltonian::operator += (const PerturbingHamiltonian& H) {
		PerturbingHamiltonian H2(H);
		for (int i = 0; i < indices.size(); i++) {//run over my indices
			int mr = indices[i].mr, mz = indices[i].mz, mphi = indices[i].mphi;
			std::vector<std::complex<double> >::iterator jt = H2.values.begin();
			for (GenFncIndices::iterator it = H2.indices.begin(); it != H2.indices.end();) {
				if (mr == (*it).mr && mz == (*it).mz && mphi == (*it).mphi) {//this index matches mine
					values[i] += (*jt);
					it = H2.indices.erase(it);
					jt = H2.values.erase(jt);
					break;//my values updated for this index
				}
				else {
					it++; jt++;//move on to the next term in H2
				}
			}
		}
		if (H2.indices.size() > 0) {//something unaccounted for
			for (int j = 0; j < H2.indices.size(); j++) {//add them into my list
				indices.push_back(H2.indices[j]);
				values.push_back(H2.values[j]);
			}
		}
		return *this;
	}

	PerturbingHamiltonian PerturbingHamiltonian::operator + (const  PerturbingHamiltonian& H) const {
		PerturbingHamiltonian H2(*this);
		H2 += H;
		return H2;
	}

#ifdef TESTIT
#include "actions_test_torus.cpp"
#endif

	coord::PosMomCyl Torus::from_toy(const Angles& thetaT) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		return from_aaT(aaT);
	}
	coord::PosMomCyl Torus::from_aaT(const ActionAngles& aaT) const {//input toy J & theta
		return TM.from_aaT(aaT);
	}
	coord::PosMomCyl Torus::from_true(const Angles& theta) const {//input true angles
		ActionAngles aaT(GF.true2toy(ActionAngles(J, theta)));//toy AAs computed from true
		return from_aaT(aaT);
	}
	//start fromT
	coord::PosCyl Torus::PosDerivJ(const Angles& thetaT,
		actions::DerivAct<coord::Cyl>& dRJ) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);

		actions::DerivAct<coord::Sph> dJ;
		actions::DerivAng<coord::Sph> dA;
		const coord::PosMomSph rtheta(Is.aa2pq(aaT, NULL, &dJ, &dA));
		double snv, csv;
		math::sincos(rtheta.theta, snv, csv);
		coord::PosUVSph uv(asinh(rtheta.r / cs.Delta), rtheta.theta, rtheta.phi, cs);
		double chu = cosh(uv.u), shu = sinh(uv.u);
		double chu2 = 2 * chu * chu - 1, shu2 = 2 * shu * chu;
		double csv2 = 2 * csv - 1, sinv2 = 2 * snv * csv;
		double dudr = 1 / (cs.Delta * chu);
		double dRdu = cs.Delta * chu * snv, dRdv = cs.Delta * shu * csv;
		double dzdu = cs.Delta * shu * csv, dzdv = -cs.Delta * chu * snv;
		//start fromT
		double n1 = (chu2 - csv2);
		double dpRdv = (2 * csv * pow_2(chu) * (-2 + csv2 + chu2)) / pow_2(n1) * rtheta.pr
			- 2 * (2 + csv2 + chu2) * snv * shu / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpzdv = -(2 + csv2 + chu2) * snv * shu2 / pow_2(n1) * rtheta.pr
			- 2 * csv * chu * (-2 + csv2 + chu2) / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpRdu = -4 * pow_2(csv) * snv * shu2 / pow_2(n1) * rtheta.pr - 2 * csv * chu *
			(-2 + csv2 + chu2) / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpzdu = -2 * csv * (csv2 * chu2 - 1) / pow_2(n1) * rtheta.pr
			+ 2 * (2 + csv2 + chu2) * snv * shu / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpRdpr = 2 * pow_2(chu) * snv / n1;
		double dpRdpt = 2 * shu * csv / (cs.Delta * n1);
		double dpzdpr = 2 * csv * chu * shu / n1;double dpzdpt = -2 * chu * snv / (cs.Delta * n1);
		//gets dw/dJT
		dRJ.dbyJr.pR = dpRdpr * dJ.dbyJr.pr + dpRdpt * dJ.dbyJr.ptheta + dpRdu * dudr * dJ.dbyJr.r
			+ dpRdv * dJ.dbyJr.theta;
		dRJ.dbyJz.pR = dpRdpr * dJ.dbyJz.pr + dpRdpt * dJ.dbyJz.ptheta + dpRdu * dudr * dJ.dbyJz.r
			+ dpRdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.pR = dpRdpr * dJ.dbyJphi.pr + dpRdpt * dJ.dbyJphi.ptheta + dpRdu * dudr * dJ.dbyJphi.r
			+ dpRdv * dJ.dbyJphi.theta;
		dRJ.dbyJr.pz = dpzdpr * dJ.dbyJr.pr + dpzdpt * dJ.dbyJr.ptheta + dpzdu * dudr * dJ.dbyJr.r
			+ dpzdv * dJ.dbyJr.theta;
		dRJ.dbyJz.pz = dpzdpr * dJ.dbyJz.pr + dpzdpt * dJ.dbyJz.ptheta + dpzdu * dudr * dJ.dbyJz.r
			+ dpzdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.pz = dpzdpr * dJ.dbyJphi.pr + dpzdpt * dJ.dbyJphi.ptheta + dpzdu * dudr * dJ.dbyJphi.r
			+ dpzdv * dJ.dbyJphi.theta;
		dRJ.dbyJr.pphi = 0.0;dRJ.dbyJz.pphi = 0.0;dRJ.dbyJphi.pphi = 1.0;
		dRJ.dbyJr.R = dRdu * dudr * dJ.dbyJr.r + dRdv * dJ.dbyJr.theta;
		dRJ.dbyJz.R = dRdu * dudr * dJ.dbyJz.r + dRdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.R = dRdu * dudr * dJ.dbyJphi.r + dRdv * dJ.dbyJphi.theta;

		dRJ.dbyJr.z = dzdu * dudr * dJ.dbyJr.r + dzdv * dJ.dbyJr.theta;
		dRJ.dbyJz.z = dzdu * dudr * dJ.dbyJz.r + dzdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.z = dzdu * dudr * dJ.dbyJphi.r + dzdv * dJ.dbyJphi.theta;
		dRJ.dbyJr.phi = dJ.dbyJr.phi;
		dRJ.dbyJz.phi = dJ.dbyJz.phi;
		dRJ.dbyJphi.phi = dJ.dbyJphi.phi;
		//dtheta_i/dthetaT_j=dJT_j/dJ_i
		math::Matrix<double> dthetadthetaT(3, 3);
		GF.dtbydtT_Jacobian(thetaT, dthetadthetaT);

		math::Matrix<double> dxdJT(3, 3);
		dxdJT(0, 0) = dRJ.dbyJr.R;dxdJT(1, 0) = dRJ.dbyJr.z;dxdJT(2, 0) = dRJ.dbyJr.phi;
		dxdJT(0, 1) = dRJ.dbyJz.R;dxdJT(1, 1) = dRJ.dbyJz.z;dxdJT(2, 1) = dRJ.dbyJz.phi;
		dxdJT(0, 2) = dRJ.dbyJphi.R;dxdJT(1, 2) = dRJ.dbyJphi.z;dxdJT(2, 2) = dRJ.dbyJphi.phi;
		//Gets dw/dJ=dw/dJT*dJt/dJ
		math::Matrix<double> Mat3(3, 3);
		math::blas_dgemm(math::CblasNoTrans, math::CblasTrans, 1.0, dxdJT, dthetadthetaT, 0.0, Mat3);
		//math::blas_dgemm(math::CblasNoTrans, math::CblasNoTrans, 1.0, dxdJT, AB, 0.0, Mat3);
		math::Matrix<double> dpdJT(3, 3);
		dpdJT(0, 0) = dRJ.dbyJr.pR;dpdJT(1, 0) = dRJ.dbyJr.pz;dpdJT(2, 0) = dRJ.dbyJr.pphi;
		dpdJT(0, 1) = dRJ.dbyJz.pR;dpdJT(1, 1) = dRJ.dbyJz.pz;dpdJT(2, 1) = dRJ.dbyJz.pphi;
		dpdJT(0, 2) = dRJ.dbyJphi.pR;dpdJT(1, 2) = dRJ.dbyJphi.pz;dpdJT(2, 2) = dRJ.dbyJphi.pphi;
		math::Matrix<double> Mat4(3, 3);
		math::blas_dgemm(math::CblasNoTrans, math::CblasTrans, 1.0, dpdJT, dthetadthetaT, 0.0, Mat4);
		dRJ.dbyJr.R = Mat3(0, 0);dRJ.dbyJr.z = Mat3(1, 0);dRJ.dbyJr.phi = Mat3(2, 0);
		dRJ.dbyJz.R = Mat3(0, 1);dRJ.dbyJz.z = Mat3(1, 1);dRJ.dbyJz.phi = Mat3(2, 1);
		dRJ.dbyJphi.R = Mat3(0, 2);dRJ.dbyJphi.z = Mat3(1, 2);dRJ.dbyJphi.phi = Mat3(2, 2);
		dRJ.dbyJr.pR = Mat4(0, 0);dRJ.dbyJr.pz = Mat4(1, 0);dRJ.dbyJr.pphi = Mat4(2, 0);
		dRJ.dbyJz.pR = Mat4(0, 1);dRJ.dbyJz.pz = Mat4(1, 1);dRJ.dbyJz.pphi = Mat4(2, 1);
		dRJ.dbyJphi.pR = Mat4(0, 2);dRJ.dbyJphi.pz = Mat4(1, 2);dRJ.dbyJphi.pphi = Mat4(2, 2);
		
		return coord::toPosCyl(uv);
	}
	//end fromT
		/* Position from toy angle plus dR/dthetaT and dpR/dthetaT at fixed J (which causes JT
		 * to vary with thetaT)
		*/
	coord::PosCyl Torus::PosDerivs(const Angles& thetaT,
		actions::DerivAngCyl& dRA, double* det) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		actions::DerivAct<coord::Sph> dJ;
		actions::DerivAng<coord::Sph> dA;
		const coord::PosMomSph rtheta(Is.aa2pq(aaT, NULL, &dJ, &dA));
		double snv, csv;
		math::sincos(rtheta.theta, snv, csv);
		coord::PosUVSph uv(asinh(rtheta.r / cs.Delta), rtheta.theta, rtheta.phi, cs);
		double chu = cosh(uv.u), shu = sinh(uv.u);
		double chu2 = 2 * chu * chu - 1, shu2 = 2 * shu * chu;
		double csv2 = 2 * csv*csv - 1, sinv2 = 2 * snv * csv;
		double dudr = 1. / (cs.Delta * chu);
		double dRdu = cs.Delta * chu * snv, dRdv = cs.Delta * shu * csv;
		double dzdu = cs.Delta * shu * csv, dzdv = -cs.Delta * chu * snv;
		//start fromT
		double n1 = (chu2-csv2);
		double dpRdv = (2 * csv * pow_2(chu) * (-2 + csv2 + chu2)) / pow_2(n1) * rtheta.pr
			- 2 * (2 + csv2 + chu2) * snv * shu / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpzdv = -(2 + csv2 + chu2) * snv * shu2 / pow_2(n1) * rtheta.pr
			- 2 * csv * chu * (-2 + csv2 + chu2) / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpRdu = -4 * pow_2(csv) * snv * shu2 / pow_2(n1) * rtheta.pr - 2 * csv * chu * 
			(-2 + csv2+ chu2) / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpzdu = -2 * csv * (csv2 * chu2 - 1) / pow_2(n1) * rtheta.pr
			+ 2 * (2 + csv2 + chu2) * snv * shu / (cs.Delta * pow_2(n1)) * rtheta.ptheta;
		double dpRdpr = 2 * pow_2(chu) * snv / n1;
		double dpRdpt = 2 * shu* csv / (cs.Delta * n1);
		double dpzdpr = 2 * csv * chu * shu / n1;double dpzdpt = -2 * chu * snv / (cs.Delta * n1);
		dRA.dbythetar.pR = dpRdpr * dA.dbythetar.pr + dpRdpt * dA.dbythetar.ptheta + dpRdu * dudr * dA.dbythetar.r
			+ dpRdv * dA.dbythetar.theta;
		dRA.dbythetaz.pR = dpRdpr * dA.dbythetaz.pr + dpRdpt * dA.dbythetaz.ptheta + dpRdu * dudr * dA.dbythetaz.r
			+ dpRdv * dA.dbythetaz.theta;
		dRA.dbythetaphi.pR = dpRdpr * dA.dbythetaphi.pr + dpRdpt * dA.dbythetaphi.ptheta + dpRdu * dudr * dA.dbythetaphi.r
			+ dpRdv * dA.dbythetaphi.theta;
		dRA.dbythetar.pz = dpzdpr * dA.dbythetar.pr + dpzdpt * dA.dbythetar.ptheta + dpzdu * dudr * dA.dbythetar.r
			+ dpzdv * dA.dbythetar.theta;
		dRA.dbythetaz.pz = dpzdpr * dA.dbythetaz.pr + dpzdpt * dA.dbythetaz.ptheta + dpzdu * dudr * dA.dbythetaz.r
			+ dpzdv * dA.dbythetaz.theta;
		dRA.dbythetaphi.pz = dpzdpr * dA.dbythetaphi.pr + dpzdpt * dA.dbythetaphi.ptheta + dpzdu * dudr * dA.dbythetaphi.r
			+ dpzdv * dA.dbythetaphi.theta;
		//end fromT
		dRA.dbythetar.R = dRdu * dudr * dA.dbythetar.r + dRdv * dA.dbythetar.theta;
		dRA.dbythetaz.R = dRdu * dudr * dA.dbythetaz.r + dRdv * dA.dbythetaz.theta;
		dRA.dbythetaphi.R = dRdu * dudr * dA.dbythetaphi.r + dRdv * dA.dbythetaphi.theta;
		dRA.dbythetar.z = dzdu * dudr * dA.dbythetar.r + dzdv * dA.dbythetar.theta;
		dRA.dbythetaz.z = dzdu * dudr * dA.dbythetaz.r + dzdv * dA.dbythetaz.theta;
		dRA.dbythetaphi.z = dzdu * dudr * dA.dbythetaphi.r + dzdv * dA.dbythetaphi.theta;
		dRA.dbythetar.phi = dA.dbythetar.phi;
		dRA.dbythetaz.phi = dA.dbythetaz.phi;
		dRA.dbythetaphi.phi = 1;

		actions::DerivAct<coord::Cyl> dRJ;
		dRJ.dbyJr.R = dRdu * dudr * dJ.dbyJr.r + dRdv * dJ.dbyJr.theta;
		dRJ.dbyJz.R = dRdu * dudr * dJ.dbyJz.r + dRdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.R = dRdu * dudr * dJ.dbyJphi.r + dRdv * dJ.dbyJphi.theta;
		dRJ.dbyJr.z = dzdu * dudr * dJ.dbyJr.r + dzdv * dJ.dbyJr.theta;
		dRJ.dbyJz.z = dzdu * dudr * dJ.dbyJz.r + dzdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.z = dzdu * dudr * dJ.dbyJphi.r + dzdv * dJ.dbyJphi.theta;
		dRJ.dbyJr.phi = dJ.dbyJr.phi;
		dRJ.dbyJz.phi = dJ.dbyJz.phi;
		dRJ.dbyJphi.phi = dJ.dbyJphi.phi;
		//start fromT
		dRJ.dbyJr.pR = dpRdu * dudr * dJ.dbyJr.r + dpRdv * dJ.dbyJr.theta + dpRdpr * dJ.dbyJr.pr + dpRdpt * dJ.dbyJr.ptheta;
		dRJ.dbyJz.pR = dpRdpr * dJ.dbyJz.pr + dpRdpt * dJ.dbyJz.ptheta + dpRdu * dudr * dJ.dbyJz.r
			+ dpRdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.pR = dpRdpr * dJ.dbyJphi.pr + dpRdpt * dJ.dbyJphi.ptheta + dpRdu * dudr * dJ.dbyJphi.r
			+ dpRdv * dA.dbythetaphi.theta;
		dRJ.dbyJr.pz = dpzdpr * dJ.dbyJr.pr + dpzdpt * dJ.dbyJr.ptheta + dpzdu * dudr * dJ.dbyJr.r
			+ dpzdv * dJ.dbyJr.theta;
		dRJ.dbyJz.pz = dpzdpr * dJ.dbyJz.pr + dpzdpt * dJ.dbyJz.ptheta + dpzdu * dudr * dJ.dbyJz.r
			+ dpzdv * dJ.dbyJz.theta;
		dRJ.dbyJphi.pz = dpzdpr * dJ.dbyJphi.pr + dpzdpt * dJ.dbyJphi.ptheta + dpzdu * dudr * dJ.dbyJphi.r
			+ dpzdv * dJ.dbyJphi.theta;
		//end fromT
		actions::DerivAng<coord::Cyl> dJA = GF.dJdt(thetaT); //dJT/dthetaT
		dRA.dbythetar.R += dRJ.dbyJr.R * dJA.dbythetar.R + dRJ.dbyJz.R * dJA.dbythetar.z;
		dRA.dbythetaz.R += dRJ.dbyJr.R * dJA.dbythetaz.R + dRJ.dbyJz.R * dJA.dbythetaz.z;
		dRA.dbythetar.z += dRJ.dbyJr.z * dJA.dbythetar.R + dRJ.dbyJz.z * dJA.dbythetar.z;
		dRA.dbythetaz.z += dRJ.dbyJr.z * dJA.dbythetaz.R + dRJ.dbyJz.z * dJA.dbythetaz.z;
		dRA.dbythetar.phi += dRJ.dbyJr.phi * dJA.dbythetar.R + dRJ.dbyJz.phi * dJA.dbythetar.z;
		dRA.dbythetaz.phi += dRJ.dbyJr.phi * dJA.dbythetaz.R + dRJ.dbyJz.phi * dJA.dbythetaz.z;
		//start fromT
		dRA.dbythetar.pR += dRJ.dbyJr.pR * dJA.dbythetar.R + dRJ.dbyJz.pR * dJA.dbythetar.z;
		dRA.dbythetaz.pR += dRJ.dbyJr.pR * dJA.dbythetaz.R + dRJ.dbyJz.pR * dJA.dbythetaz.z;
		dRA.dbythetar.pz += dRJ.dbyJr.pz * dJA.dbythetar.R + dRJ.dbyJz.pz * dJA.dbythetar.z;
		dRA.dbythetaz.pz += dRJ.dbyJr.pz * dJA.dbythetaz.R + dRJ.dbyJz.pz * dJA.dbythetaz.z;
		dRA.dbythetar.pphi = dJA.dbythetar.pphi;dRA.dbythetaz.pphi = dJA.dbythetaz.pphi;dRA.dbythetaphi.pphi = dJA.dbythetaphi.pphi;
		//end fromT
		if (det) {//assume only index.mphi=0 non-zero
			(*det) = dRA.dbythetar.R * dRA.dbythetaz.z - dRA.dbythetar.z * dRA.dbythetaz.R;
		}
		return coord::toPosCyl(uv);
	}
	void Torus::zSoS(std::vector<double>& Rs, std::vector<double>& vRs, const int N,
		double& Rmin, double& Rmax, double& Vmax) const {
		const double tol = 1e-5;
		Rmin = 1e10, Rmax = 0, Vmax = 0;
		for (int i = 0; i < N; i++) {
			double thetaT_r = 2 * M_PI / (double)N * (-N / 2 + i);
			CrossingFinder CF(this, thetaT_r);
			double thetaT_z, dtheta = .1, th_min = -.5 * M_PI,
				th_max = th_min + dtheta, z_min, z_max;
			CF.evalDeriv(th_min, &z_min); CF.evalDeriv(th_max, &z_max);
			while (z_min * z_max > 0) {// plod round looking for sign change
				th_min = th_max; z_min = z_max; th_max += dtheta;
				CF.evalDeriv(th_max, &z_max);
			}
			thetaT_z = math::findRoot(CF, th_min, th_max, tol);
			coord::PosMomCyl Rz(from_toy(Angles(thetaT_r, thetaT_z, 0)));
			if (Rz.pz < -.0001) {// keep going round
				do {
					th_min = th_max; z_min = z_max;
					th_max += dtheta; CF.evalDeriv(th_max, &z_max);
				} while (z_min * z_max > 0);
				thetaT_z = math::findRoot(CF, th_min, th_max, tol) - 2 * M_PI;
				Rz = coord::PosMomCyl(from_toy(Angles(thetaT_r, thetaT_z, 0)));
			}
			Rs.push_back(Rz.R); vRs.push_back(Rz.pR);
			Rmin = fmin(Rmin, Rz.R); Rmax = fmax(Rmax, Rz.R); Vmax = fmax(Vmax, fabs(Rz.pR));
		}
		Rs.push_back(Rs[0]); vRs.push_back(vRs[0]);
	}
	std::vector<std::pair<coord::PosVelCyl, double> > Torus::orbit(const Angles& theta0, double dt, double T) const {
		std::vector<std::pair<coord::PosVelCyl, double> > traj;
		double t = 0;
		while (t < T) {
			actions::Angles theta(theta0.thetar + freqs.Omegar * t, theta0.thetaz + freqs.Omegaz * t,
				theta0.thetaphi + freqs.Omegaphi * t);
			traj.push_back(std::pair<coord::PosVelCyl, double>(coord::toPosVelCyl(from_true(theta)), t));
			t += dt;
		}
		return traj;
	}

	// Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise. If the 
	// torus passes through the point given, this happens four times, in each case
	// with a different velocity, but only two of these are independent:
	// theta_r -> -theta_r with theta_z -> Pi-theta_z leaves (R,z) and J^T fixed
	// but changes the sign of both velocities. |d(x,y,z)/d(theta_r,theta_z,theta_phi)|
	// is returned. The latter vanishes on the edge of the
	// orbit, such that its inverse, the density of the orbit, diverges there
	// (that's the reason why the density itself is not returned).

	bool Torus::containsPoint(const coord::PosCyl& p, std::vector<Angles>& As,
		std::vector<coord::VelCyl>& Vs,
		std::vector<double>& Jacobs,
		std::vector<DerivAngCyl>* dRdtheta,
		const double tol) const {
		coord::PosMomCyl peri(from_true(Angles(0, .5 * M_PI, 0))), apo(from_true(Angles(M_PI, 0, 0))),
			top(from_true(Angles(M_PI, .5 * M_PI, 0)));
		double Rmin = .95 * peri.R, Rmax = 1.05 * apo.R, zmax = 1.05 * fabs(top.z);
		if (p.R<Rmin || p.R>Rmax || fabs(p.z) > zmax) return false;
		locFinder LF(*this, p);
		double tolerance = 1e-8;
		double params[2] = { 1,1 }, result[2], dist, det;
		int maxNumIter = 200;
		coord::PosMomCyl P1; Angles A1, Atrue;
		DerivAngCyl dA;
		int done, kmax = 30, nfail = 0;

		while (As.size() < 4 && nfail < kmax) {
			double kount = 0;
			for (int k = 0; k < kmax; k++) {
				done = math::nonlinearMultiFit(LF, params, tolerance, maxNumIter, result);

				A1 = Angles(math::wrapAngle(result[0]), math::wrapAngle(result[1]), 0.0);
				P1 = from_toy(A1); P1.phi = p.phi;
				//from_toy set phi=0
				A1.thetaphi = TM.pq2aa(P1).thetaphi;

				ActionAngles aaT = TM.pq2aa(P1);
				dist = sqrt(pow_2(p.R - P1.R) + pow_2(p.z - P1.z));
				if (dist < 2 * tol) {
					Atrue = GF.trueA(A1);
					if (is_new(Atrue, As)) {
						break;
					}
				} //else	printf("Max in containsPoint: %d %g\n", done, dist);
				params[0] += .3; params[0] = math::wrapAngle(params[0]);
				params[1] += .7;   params[1] = math::wrapAngle(params[1]);
				kount++;
				if (kount == kmax && As.size() == 0) return false;
			}
			if (kount < kmax) {
				math::Matrix<double> M(3, 3);//to hold dtheta/dthetaT
				As.push_back(Atrue); PosDerivs(A1, dA, &det);
				Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi / P1.R));
				Jacobs.push_back(fabs(det / GF.dtbydtT_Jacobian(A1, M)));
				if (dRdtheta) {
					assemble(dA, M); dRdtheta->push_back(dA);
				}
				A1.thetar = -A1.thetar; A1.thetaz = M_PI - A1.thetaz;
				P1 = from_toy(A1); P1.phi = p.phi;
				A1.thetaphi = TM.pq2aa(P1).thetaphi;
				Atrue = GF.trueA(A1);
				As.push_back(Atrue); PosDerivs(A1, dA, &det);
				Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi / P1.R));
				Jacobs.push_back(fabs(det / GF.dtbydtT_Jacobian(A1, M)));
				if (dRdtheta) {
					assemble(dA, M); dRdtheta->push_back(dA);
				}
			}
			else nfail++;
		}
		if (nfail >= kmax) printf("containsPoint error at Rz (%f %f) - %d angles \n",
			p.R, p.z, As.size());
		return true;
	}
	double Torus::density(const coord::PosCyl& Rz) const {
		std::vector<Angles> As; std::vector<coord::VelCyl> Vs;
		std::vector<double> Jacobs;
		//const double tol = 1e-6;
		if (!containsPoint(Rz, As, Vs, Jacobs)) return 0;
		double rho = 0;
		for (int i = 0; i < As.size(); i++)
			rho += 1 / Jacobs[i];
		return rho;
	}
	void Torus::write(FILE* ofile) const {
		fprintf(ofile, "%g %g %g %g %g %g %g %g %g %g\n",
			J.Jr, J.Jz, J.Jphi, freqs.Omegar, freqs.Omegaz, freqs.Omegaphi,
			E, cs.Delta, Is.Js, Is.b);
		GF.write(ofile);
	}
	void Torus::read(FILE* ifile) {
		double Delta, Js, b;
		fscanf_s(ifile, "%g %g %g %g %g %g %g %g %g %g\n",
			&J.Jr, &J.Jz, &J.Jphi, &freqs.Omegar, &freqs.Omegaz, &freqs.Omegaphi,
			&E, &Delta, &Js, &b);
		cs = coord::UVSph(Delta);
		Is = Isochrone(Js, b);
		GF.read(ifile);
	}
	/* Operators on Tori
	 */
	Torus& Torus::operator *= (const double a) {
		J *= a; freqs *= a; GF *= a; Is *= a; cs = coord::UVSph(cs.Delta * a); E* a;
		TM = ToyMap(cs, Is);
		return *this;
	}
	Torus& Torus::operator += (const Torus& T) {
		J += T.J; freqs += T.freqs; GF += T.GF; Is += T.Is;
		cs = coord::UVSph(cs.Delta += T.cs.Delta); E + T.E;
		TM = ToyMap(cs, Is);
		return *this;
	}
	const Torus Torus::operator * (const double a) const {
		Torus T2(J * a, freqs * a, GF * a, Is * a, coord::UVSph(cs.Delta * a), E * a);
		return T2;
	}
	const Torus Torus::operator + (const Torus& T) const {
		Torus T2(J + T.J, freqs + T.freqs, GF + T.GF, Is + T.Is, coord::UVSph(cs.Delta + T.cs.Delta), E + T.E);
		return T2;
	}
	/* Same foe eTori
	 */
	eTorus& eTorus::operator *= (const double a) {
		J *= a; freqs *= a; GF *= a; Is *= a; cs = coord::UVSph(cs.Delta * a); E* a; pH* a;
		TM = ToyMap(cs, Is);
		return *this;
	}
	eTorus& eTorus::operator += (const eTorus& T) {
		J += T.J; freqs += T.freqs; GF += T.GF; Is += T.Is;
		cs = coord::UVSph(cs.Delta += T.cs.Delta); E + T.E; pH + T.pH;
		TM = ToyMap(cs, Is);
		return *this;
	}
	const eTorus eTorus::operator * (const double a) const {
		eTorus T2(J * a, freqs * a, GF * a, Is * a, coord::UVSph(cs.Delta * a), E * a, pH * a);
		return T2;
	}
	const eTorus eTorus::operator + (const eTorus& T) const {
		eTorus T2(J + T.J, freqs + T.freqs, GF + T.GF, Is + T.Is, coord::UVSph(cs.Delta + T.cs.Delta),
			E + T.E, pH + T.pH);
		return T2;
	}
	/* Gather all terms in the pH that are harmonics of the specified line
	 */
	std::vector<std::complex<double> > PerturbingHamiltonian::get_hn(const GenFncIndex& I,
		std::vector<float>& multiples) const {
		std::vector<std::complex<double> > Hs;
		printf("Looking for (%d %d %d)\n", I.mr, I.mz, I.mphi);
		for (int i = 0; i < indices.size(); i++) {
			GenFncIndex In(indices[i]);
			if ((I.mr == 0 && In.mr != 0) || (I.mr != 0 && In.mr == 0)) continue;
			if ((I.mz == 0 && In.mz != 0) || (I.mz != 0 && In.mz == 0)) continue;
			if ((I.mphi == 0 && In.mphi != 0) || (I.mphi != 0 && In.mphi == 0)) continue;
			//now either matching zero indices or both non-zero
			std::vector<double> Rats;
			if (I.mr != 0) Rats.push_back((double)I.mr / (double)In.mr);
			if (I.mz != 0) Rats.push_back((double)I.mz / (double)In.mz);
			if (I.mphi != 0) Rats.push_back((double)I.mphi / (double)In.mphi);
			if (Rats[1] != Rats[0]) continue;
			if ((Rats.size() == 3) && (Rats[2] != Rats[0])) continue;
			printf("%3d %2d %2d %g %f at rank %d\n", In.mr, In.mz, In.mphi, math::modulus(values[i]),
				math::arg(values[i]) / M_PI, i);
			Hs.push_back(values[i]);
			if (Rats[0] >= 1) multiples.push_back(Rats[0]);
			else multiples.push_back(1 / Rats[0]);
		}
		return Hs;
	}
	TMfitter::TMfitter(const potential::BasePotential& pot,
		std::vector<std::pair<coord::PosMomCyl, double> >& _traj,
		double _pphi) : traj(_traj), pphi(_pphi) {
		xmin = 1e6; ymax = 0;
		double ymin = 1e6, xmax = 0, pxmin = 1e6, pxmax = 0, pymin = 1e6, pymax = 0;
		for (int i = 1; i < traj.size(); i++) {//Find where crosses axes
			if (traj[i].first.R * traj[i - 1].first.R <= 0) {//x axis
				double f = traj[i].first.R / (traj[i].first.R - traj[i - 1].first.R);
				double y = (1 - f) * traj[i].first.z + f * traj[i - 1].first.z;
				double px = (1 - f) * traj[i].first.pR + f * traj[i - 1].first.pR;
				ymin = fmin(ymin, fabs(y)); ymax = fmax(ymax, fabs(y));
				pxmin = fmin(pxmin, fabs(px)); pxmax = fmax(pxmax, fabs(px));
			}
			if (traj[i].first.z * traj[i - 1].first.z <= 0) {//z axis
				double f = traj[i].first.z / (traj[i].first.z - traj[i - 1].first.z);
				double x = (1 - f) * traj[i].first.R + f * traj[i - 1].first.R;
				double py = (1 - f) * traj[i].first.pz + f * traj[i - 1].first.pz;
				xmin = fmin(xmin, fabs(x)); xmax = fmax(xmax, fabs(x));
				pymin = fmin(pymin, fabs(py)); pymax = fmax(pymax, fabs(py));
			}
		}
		//xbar, ybar estimated axes of underlying loop orbit
		xbar = .5 * (xmin + xmax); double ybar = .5 * (ymin + ymax);
		Delta2 = (ybar * ybar - xbar * xbar);
		double Phi; coord::GradCyl grad;
		pot.eval(coord::PosCyl(0, ymax, 0), &Phi, &grad);
		Frat = grad.dz;
		pot.eval(coord::PosCyl(xmin, 0, 0), &Phi, &grad);
		Frat /= grad.dR;
		double pxbar = pphi > 0 ? -.5 * (pxmin + pxmax) : .5 * (pxmin + pxmax),
			pybar = pphi > 0 ? .5 * (pymin + pymax) : -.5 * (pymin + pymax);
		aPT = 0.25 * pphi * (1 / (ybar * pybar) + 1 / (xbar * pxbar));
		bPT = 0.125 * pphi * (1 / (ybar * pybar) - 1 / (xbar * pxbar)) - 0.25;
		printf("xbar %f ybar %f ", xbar, ybar);
		printf("pxbar %f pybar %f\n", pxbar, pybar);
		printf("Delta2 %f Frat %f a %f b %f\n", Delta2, Frat, aPT, bPT);
	}
	//Function with root where we have the right force ratio
	double TMfitter::value(double Js) const {
		double g2 = pow((fabs(pphi) + sqrt(pphi * pphi + 4 * Js * Js)) / (2 * Js), 4) - 1;
		double bIso = xbar / sqrt(g2);
		double ap = sqrt(ymax * ymax + bIso * bIso), am = sqrt(xmin * xmin + bIso * bIso);
		//	printf("Js pphi g2 b %f %g %g %g\n",Js,pphi,g2,b);
		return Frat - (ymax * am * pow_2(bIso + am)) / (xmin * ap * pow_2(bIso + ap));
	}
	//Solve for Js and bIso
	std::vector<double> TMfitter::fitTM() const {
		const double Jsmin = .01, Jsmax = 10;
		double Js = math::findRoot(*this, Jsmin, Jsmax, 1e-5);
		double g2 = (pow((fabs(pphi) + sqrt(pphi * pphi + 4 * Js * Js)) / (2 * Js), 4) - 1);
		double bIso = xbar / sqrt(g2);
		std::vector<double> ans(5);
		ans[0] = sqrt(Delta2); ans[1] = aPT; ans[2] = bPT;
		ans[3] = Js; ans[4] = bIso;
		return ans;
	}

	TorusGenerator::TorusGenerator(const potential::BasePotential& _pot, const double _tol) :
		pot(_pot), defaultTol(_tol), invPhi0(1. / _pot.value(coord::PosCyl(0, 0, 0))), tmax(250) {
		std::vector<double> gridR = potential::createInterpolationGrid(pot, ACCURACY_INTERP2);
		int sizeL = gridR.size(), sizeXi = 20;
		printf("Preparing TorusGenerator...");
		math::ScalingSemiInf Lscale;
		std::vector<double> gridL(sizeL);
		std::vector<double> gridLscaled(sizeL);
		std::vector<double> gridXi(sizeXi);
		for (int i = 0; i < sizeL; i++) {
			gridL[i] = gridR[i] * potential::v_circ(pot, gridR[i]);
			//		gridLscaled[i]=i/(double)(sizeL-1);
			//		gridL[i]=math::unscale(Lscale,gridLscaled[i]);
		}
		for (int i = 0; i < sizeXi; i++)
			gridXi[i] = i / (double)(sizeXi - 1);//EV wld call this Xiscaled
		math::Matrix<double> grid2dD(sizeL, sizeXi);
		math::Matrix<double> grid2dR(sizeL, sizeXi);
		createGridFocalDistance(pot, gridL, gridXi, grid2dD, grid2dR);
		interpD = math::LinearInterpolator2d(gridL, gridXi, grid2dD);
		interpR = math::LinearInterpolator2d(gridL, gridXi, grid2dR);
		printf("done\n");
	}
	double TorusGenerator::Hamilton(const Torus& T, const potential::BasePotential* ePot, const Angles& theta)
	{
		coord::PosMomCyl Rz(T.from_true(theta));
		double H = .5 * (pow_2(Rz.pR) + pow_2(Rz.pz) + pow_2(Rz.pphi / Rz.R)) + pot.value(Rz);
		coord::PosCyl pos(Rz.R, Rz.z, Rz.phi);
		return ePot ? H + ePot->value(pos) : H;
	}
	PerturbingHamiltonian TorusGenerator::get_pH(const Torus& T, int nf, bool ifp,
		const potential::BasePotential* ePot) {//Fourier analyses H
		int nfr = nf, nfz = nf;
		int nfp = ePot ? nf / 4 : 1;
		double N = (nfp * nfr * nfz);
		double* h = new double[nfp * nfr * nfz];
		Angles thetas;
		double dtr = 2 * M_PI / (double)nfr, dtz = 2 * M_PI / (double)nfz, dtp = M_PI / (double)nfp;
#pragma omp parallel for schedule(dynamic) 
		for (int k = 0; k < nfp; k++) {
			thetas.thetaphi = k * dtp;//if !ePot sticks at 0
			for (int i = 0; i < nfr; i++) {
				int i1 = nfr - i;
				thetas.thetar = i * dtr;
				for (int j = 0; j < nfz; j++) {
					int j1 = nfz - j, j2 = nfz / 2 - j, j3 = nfz / 2 + j;
					thetas.thetaz = j * dtz;
					double a = Hamilton(T, ePot, thetas);
					h[nfr * nfz * k + nfz * i + j] = a;
					/*
					if(j3<nfz) h[nfr*nfz*k+i*nfz+j3]=a;//N-S symmetry
					if(i1!=i && i1<nfr){
						if(j1<nfz) h[nfr*nfz*k+i1*nfz+j1]=a;//time-reverse symmetry
						if(j2<nfz) h[nfr*nfz*k+i1*nfz+j2]=a;//both symmetries
					}*/
				}
			}
		}
		//need speq to return cpts at Nyquist vals m_phi, which we won't use 
		math::Matrix<double> speq(nfr, 2 * nfz);
		rlft3(h, speq, nfp, nfr, nfz, 1);
		std::vector<double> Hmods;
		GenFncIndices Hindices;
		int ntop = 0;
		double hmax = 0;
		for (int k = 0; k < nfp; k++) {//find largest perturbing terms
			for (int i = 0; i < nfr; i++) {
				for (int j = 0; j < nfz / 2; j++) {
					double s = sqrt(pow_2(h[nfr * nfz * k + nfz * i + 2 * j])
						+ pow_2(h[nfr * nfz * k + nfz * i + 2 * j + 1])) / N;
					if (ntop<tmax || s>hmax) {
						hmax = insertLine(ntop, tmax, s, GenFncIndex(i, j, k),
							Hmods, Hindices);
					}
				}
			}
		}
		ifp = false;
		std::vector<std::complex<double> > Hvalues;
		for (int i = 0; i < ntop; i++) {
			Hvalues.push_back(std::complex<double>(h[nfr * nfz * Hindices[i].mphi + nfz * Hindices[i].mr + 2 * Hindices[i].mz] / N,
				h[nfr * nfz * Hindices[i].mphi + nfz * Hindices[i].mr + 2 * Hindices[i].mz + 1] / N));
		}
		if (ifp) printf("Terms in perturbing H (*100)\n");
		for (int i = 0; i < ntop; i++) {
			//		printf("%d %d %d ",Hindices[i].mr, Hindices[i].mz, Hindices[i].mphi);
			if (Hindices[i].mphi > nfp / 2) Hindices[i].mphi -= nfp;
			if (Hindices[i].mr > nfr / 2)   Hindices[i].mr -= nfr;
			//		Hindices[i].mz*=2;
			//		if(Hindices[i].mz>nfz/2)   Hindices[i].mz   -= nfz;//won't happen!
			if (ifp) printf("(%3d %3d %3d) (%g %g)\n",
				Hindices[i].mr, Hindices[i].mz, Hindices[i].mphi,
				100 * math::modulus(Hvalues[i]), math::arg(Hvalues[i]));
		}
		delete[] h;
		return PerturbingHamiltonian(Hindices, Hvalues);
	}

	double TorusGenerator::getRsh(Actions& J) {
		const double L = J.Jz + fabs(J.Jphi);
		return interpR.value(L, J.Jz / L);
	}
	/* setConsts fixes Jscale Delta, Rsh and Is
	 */
	void TorusGenerator::setConsts(actions::Actions J, double& Jscale, double& freqScale,
		double& Rsh, Isochrone& Is, coord::UVSph& cs) const {
		const double relToler = 1e-4;
		const double L = fabs(J.Jphi) + J.Jz, Xi = J.Jz / L;
		const double Jtot = L + J.Jr;
		Jscale = J.Jr + J.Jz;
		Iso ISO(L, J.Jr);
		double Delta = interpD.value(L, Xi);
		cs = coord::UVSph(Delta);
		Rsh = interpR.value(L, Xi);
		freqScale = potential::v_circ(pot, Rsh) / Rsh; //frequency scale set
		//For any Js b=Rsh/ISO.g(Js) f_apo_peri=ISO.f(b.Js,e) & where
		//this matches F_apo_peri in pot we pick Js and b
		double Jsmax = 1.1 * Jtot, Jsmin = .1 * Jsmax;
		JsFinder JF(pot, ISO, Rsh);
		double val1, val2;
		JF.evalDeriv(Jsmin, &val1); JF.evalDeriv(Jsmax, &val2);
		while (val1 > 0) {
			Jsmin *= .75; JF.evalDeriv(Jsmin, &val1);
		}
		while (val2 < 0) {
			Jsmax *= 1.5; JF.evalDeriv(Jsmax, &val2);
		}
		double Js_iso = math::findRoot(JF, Jsmin, Jsmax, relToler);
		//	printf("%g %g %g %g\n",val1,val2,Js_iso,ISO.g(Js_iso));
		double b_iso = Rsh / ISO.g(Js_iso);
		Is = Isochrone(Js_iso, b_iso);
#ifdef PLT
		printf("Js: %6.3f, b: %6.3f, Delta: %6.3f, Rsh: %6.3f, freqScale: %6.3f\n",
			Js_iso, b_iso, Delta, Rsh, freqScale);
		const int ns = 40;
		double F[2], Js[ns], f_apo_peri[ns], F_apo_peri[ns];
		for (int i = 0; i < ns; i++) {
			Js[i] = Jsmin + i * (Jsmax - Jsmin) / (double)(ns - 1);
			double b = Rsh / ISO.g(Js[i]);
			double e;
			f_apo_peri[i] = ISO.f(b, Js[i], e);
			double c = ISO.cob(Js[i]) * b;
			for (int k = -1; k < 2; k += 2) {
				double u = 1 + k * e;
				double Phi, r = c * sqrt(u * (u + 2 * b / c));
				if (k == -1) printf("(%f ", r); else printf("%f ", r);
				coord::PosCyl Rz(r, 0, 0); coord::GradCyl grad;
				pot.eval(Rz, &Phi, &grad);
				F[(k + 1) / 2] = grad.dR;
			}
			F_apo_peri[i] = F[1] / F[0];
			printf("%f) ", b);
		}
		mgo::plt pl;
		pl.new_plot(Jsmin, Jsmax, 0, 1, "J\\ds", "F\\da/F\\dp");
		pl.connect(Js, f_apo_peri, ns);
		pl.setcolour("red");
		pl.connect(Js, F_apo_peri, ns);
		pl.relocate(Jsmin, .1);pl.label(" true \\gF");
		pl.grend();
#endif
	}
	Torus TorusGenerator::fitTorus(const Actions& J, const double tighten) const {
		double tolerance = 1e-8;//controls optimisation of the given Sn
		double Jscale, freqScale, Rsh;
		Isochrone Is;
		coord::UVSph cs;
		setConsts(J, Jscale, freqScale, Rsh, Is, cs);
		std::vector<double> params;
		int nrmax = 7;
		int nzmax = 6;// nzmax must be even
		double Hbar, Hdisp = 1e20;
		double Hbar2 = 1e20;
		bool converged = false;
		int Loop = 0, MaxLoop = 10, maxNumIter = 10;
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		GenFncFracs fracs = {};
		GenFncFracs fracs2;
		for (int i = 0;i <= 2;i += 2) {
			fracs2.push_back(GenFncFrac(i));
		}
		//fracs2 = {};
		GenFncDerivs derivs;
		params.resize(indices.size(), 0);
		int np = params.size() - 1;

		double tol = defaultTol * tighten;
		std::vector<double> rep;
		do {
			GenFncFitSeries GFFS(indices, fracs, J);
			torusFitter TF(J, pot, Is, cs, GFFS);
			//printf("TF: %d angles, %d parameters ", TF.numValues(), TF.numVars());
			//printf("GFFS: %d %d %d\n", GFFS.numPoints(), GFFS.numParams(), GFFS.numFracs());
			try {
				int numIter = math::nonlinearMultiFit(TF, &params[0], tolerance, maxNumIter, &params[0]);
				Hdisp = sqrt(TF.computeHamiltonianDisp(params, Hbar));
				rep.push_back(Hdisp);
				converged = (Hdisp < tol * freqScale * Jscale);
			}
			catch (std::exception& e) {
				std::cout << "Exception in fitTorus: " << e.what() << '\n';
			}
			if (converged) break;
			Loop++;
			if (Loop < MaxLoop) {
				for (int k = 0;k < fracs2.size();k++){
					int mz = fracs2[k].mz;double B, disp;
					GenFnc G2(indices, params, derivs);
					double b = G2.TermRatio(4, mz, B, &disp);
					if (disp < 1e-4) {
						int maxr = 0;
						for (int k = 0;k < indices.size();k++) {
							if (indices[k].mz == mz) maxr = std::max<int>(maxr, indices[k].mr);
						}
						if (maxr > 10) {
							//printf("maxr:%d\n", maxr);
							fracs2[k].krmin = maxr + 1;
							fracs.push_back(fracs2[k]);
							params.resize(indices.size() + 2 * fracs.size(), 0);
							np = params.size() - 1;
							params[np - 1] = B*pow(b,(maxr+1)); params[np] = atanh(b);
							//remove this fracs from list
							fracs2.erase(fracs2.begin()+k);
						}
					}
				}
				GFFS=GenFncFitSeries(indices, fracs, J);
				indices = GFFS.expand(params);
				
			}
		} while (Loop < MaxLoop);
		if (!converged) {
			printf("fitTorus failed to converge for target %g\n", tol * freqScale * Jscale);
			for (int i = 0; i < rep.size(); i++) printf("%g ", rep[i]); printf("\n");
		}
		//printf("Converged after %d loops with dispersion %g\n", Loop,Hdisp);
		//GenFncFitSeries GFFS2(indices, fracs2, J);
		//torusFitter TF2(J, pot, Is, cs, GFFS2);
		GenFncFitSeries GFFS(indices, fracs, J);
		torusFitter TF(J, pot, Is, cs, GFFS);
#ifdef TESTIT
		test_it(J, params);
#endif
		Frequencies freqs;
		//GenFncDerivs derivs;
		Hdisp = TF.fitAngleMap(&params[0], Hbar, freqs, derivs);
		
		//std::vector<double> params2(indices.size());
		GenFnc G(indices, params, derivs, fracs);
		return Torus(J, freqs, G, Is, cs, Hbar);
	}
	Torus TorusGenerator::fitFrom(const Actions& J, const Torus T, const double tighten) const {
		double tolerance = 1e-8;//controls optimisation of the given Sn
		double freqScale = .5 * (T.freqs.Omegar + T.freqs.Omegaphi);
		double Jscale = J.Jr + J.Jz;
		Isochrone Is = T.Is;
		coord::UVSph cs = T.cs;
		GenFncIndices indices = T.GF.indices;
		std::vector<double> params = T.GF.values;
		GenFncFracs fracs = {};
		GenFncFracs fracs2;
		GenFncDerivs derivs;
		for (int i = 0;i <= 2;i += 2) {
			fracs2.push_back(GenFncFrac(i));
		}
		params.resize(indices.size() , 0);
		int np = params.size() - 1;
		for (int j = 0; j < fracs.size(); j++) {
			params[np - 2 * j - 1] = 1e-3; params[np - 2 * j] = .5;//B, b
		}
		double Hbar, Hdisp = 1e20;
		bool converged = false;
		int Loop = 0, MaxLoop = 10, maxNumIter = 10;
		double tol = defaultTol * tighten;
		std::vector<double> rep;
		do {
			GenFncFitSeries GFF(indices,fracs, J);
			params.resize(indices.size());
			torusFitter TF(J, pot, Is, cs, GFF);
			try {
				int numIter = math::nonlinearMultiFit(TF, &params[0], tolerance, maxNumIter, &params[0]);
				Hdisp = sqrt(TF.computeHamiltonianDisp(params, Hbar));
				rep.push_back(Hdisp);
				converged = (Hdisp < tol * freqScale * Jscale);
				if (TF.giveNANfrac() != 0)
					printf("Fraction %7.3f of H vals NANs\n", TF.giveNANfrac());
			}
			catch (std::exception& e) {
				std::cout << "Exception in fitTorus: " << e.what() << '\n';
			}
			if (converged) break;
			Loop++;
			if (Loop < MaxLoop) {
				indices = GFF.expand(params);
				for (int k = 0;k < fracs2.size();k++) {
					int mz = fracs2[k].mz;double B, disp;
					GenFnc G2(indices, params, derivs);
					double b = G2.TermRatio(4, mz, B, &disp);
					if (disp < 1e-4) {
						int maxr = 0;
						for (int k = 0;k < indices.size();k++) {
							if (indices[k].mz == mz) maxr = std::max<int>(maxr, indices[k].mr);
						}
						if (maxr > 10) {
							printf("maxr:%d\n", maxr);
							fracs2[k].krmin = maxr + 1;
							fracs.push_back(fracs2[k]);
							params.resize(indices.size() + 2 * fracs.size(), 0);
							np = params.size() - 1;
							params[np - 1] = B * pow(b, (maxr + 1)); params[np] = atanh(b);
							//remove this fracs from list
							fracs2.erase(fracs2.begin() + k);
						}
					}
				}
			}
		} while (Loop < MaxLoop);
		if (!converged) {
			printf("fitTorus failed to converge for target %g\n", tol * freqScale * Jscale);
			for (int i = 0; i < rep.size(); i++) printf("%g ", rep[i]); printf("\n");
		}
		GenFncFitSeries GFF(indices,fracs, J);
		torusFitter TF(J, pot, Is, cs, GFF);
#ifdef TESTIT
		test_it(J, params);
#endif
		Frequencies freqs;
		Hdisp = TF.fitAngleMap(&params[0], Hbar, freqs, derivs);
		//printf("Hbar: %f, freqs: %f\n", Hbar, freqs.Omegaz);
		GenFnc G(indices, params, derivs,fracs);
		return Torus(J, freqs, G, Is, cs, Hbar);
	}
	//*/
	eTorus TorusGenerator::fiteTorus(const Actions& J, const double tol, const potential::BasePotential* ePhi) {
		Torus T = fitTorus(J, tol);
		int nf = 128;
		PerturbingHamiltonian pH(get_pH(T, nf, true, ePhi));
		return eTorus(T, pH);
	}
	eTorus TorusGenerator::fiteTorus(const Actions& J, const potential::BasePotential* ePhi) {
		Torus T = fitTorus(J);
		int nf = 128;
		PerturbingHamiltonian pH(get_pH(T, nf, true, ePhi));
		return eTorus(T, pH);
	}
	double TorusGenerator::getDelta(double& L, double& Xi) {
		return interpD.value(L, Xi);
	}
	double TorusGenerator::getDelta(actions::Actions& J) {
		double L = fabs(J.Jphi) + J.Jz, Xi = J.Jz / L;
		return interpD.value(L, Xi);
	}
	std::vector<Torus> TorusGenerator::constE(const double Jrmin, const Actions& Jstart, const int Nsteps) {
		double fac = exp(-log(Jstart.Jr / Jrmin) / (Nsteps - 1));
		Actions Jnext(Jstart);
		std::vector<Torus> Tgrd;
		Torus T(fitTorus(Jnext));
		Tgrd.push_back(T);
		for (int i = 1; i < Nsteps; i++) {
			double frat = T.freqs.Omegar / T.freqs.Omegaz;
			double dJr = (1 - fac) * Jnext.Jr;
			Jnext.Jr -= dJr; Jnext.Jz += dJr * frat;
			//printf("Next Jr: %f\n",Jnext.Jr);
			Torus T1(fitTorus(Jnext));
			T1 *= 0.5; T1 += T * 0.5;
			double frat1 = T1.freqs.Omegar / T1.freqs.Omegaz;
			Jnext.Jz = Jnext.Jz + dJr * frat - dJr * frat1;
			T = fitTorus(Jnext);
			Tgrd.push_back(T);
		}
		return Tgrd;
	}
	EXP ActionAngles ActionFinderTG::actionAngles2(const coord::PosVelCyl& point, Frequencies* freq) const {
		double E = potential::totalEnergy(*pot, point);
		if (E > 0) {
			printf("Energy is: %f and positive so no orbit\n", E);
			return ActionAngles(actions::Actions(NAN, NAN, NAN), actions::Angles(NAN, NAN, NAN));
		}
		const double tol = .5e-5;
		double phi0 = point.phi, last_diff = 1e6;
		while (fabs(phi0) > M_PI) phi0 += phi0 > M_PI ? -M_PI : M_PI;
		coord::PosMomCyl P0(point.R, point.z, phi0, point.vR, point.vz, point.vphi * point.R);
		Angles trueA;
		Actions J(AF.actions(point));
		Torus T(TG.fitTorus(J));
		double b = T.Is.b;
		double del = T.cs.Delta;
		double r2 = point.R * point.R + point.z * point.z;
		double r = sqrt(0.5 * ((r2 + del * del) + sqrt(pow_2(r2 + del * del) + 4 * del * del * point.R * point.R)));
		coord::GradCyl dR;
		coord::PosCyl x1 = coord::PosCyl(point);
		(*pot).eval(x1, NULL, &dR);
		double Js = sqrt((dR.dR * x1.R / r + dR.dz * r * x1.z / (r * r + del * del)) * b * pow_2(b + sqrt(r * r + b * b)) * sqrt(r * r + b * b) / r);
		
		//std::cout << Js << " " << T.Is.Js << '\n';
		//T.Is.Js = Js;
		//T = TG.fitFrom(J, T);
		int kount = 0;
		ActionAngles aaT = T.TM.pq2aa(P0);
		Angles tT(aaT);
		Actions JT(aaT);
		Actions JT2 = T.GF.toyJ(J, tT);
		
		double diff = 1e3;
		while (kount < 10) {
			if (kount==1) {
				T = Torus(TG.fitTorus(J));
				aaT= T.TM.pq2aa(P0);
				JT = Actions(aaT);
				tT = Angles(aaT);
			}
			else if (kount > 0) T = TG.fitFrom(J, T);
			Actions JT1 = T.GF.toyJ(J, tT);
			std::vector<double> df = {JT.Jr-JT1.Jr,JT.Jz - JT1.Jz,0.0};
			
			math::Matrix<double> dthetadthetaT(3, 3);
			T.GF.dtbydtT_Jacobian(tT, dthetadthetaT);
			math::Matrix<double> Mat3(3, 3);
			math::LUDecomp LUM(dthetadthetaT);
			math::Matrix<double> inv=LUM.inverse(3);
			std::vector<double> dJt(3,0.0);
			math::blas_dgemv(math::CblasTrans, 1.0, inv, df, 0.0, dJt);
			//std::vector<double> dJt = LUM.solve(df);
			diff = sqrt(pow_2(dJt[0]) + pow_2(dJt[1]));
			J.Jr += dJt[0]; J.Jz += dJt[1];J.Jphi += dJt[2];
			if (sqrt(pow_2(dJt[0]) + pow_2(dJt[1])) < tol) {
				if(kount>0) {
					if (freq) *freq = T.freqs;
					trueA = T.GF.trueA(tT);
					break;
				}
			}
			kount++;
		}
		return ActionAngles(J, trueA);
	}
	EXP ActionAngles ActionFinderTG::actionAngles(const coord::PosVelCyl& point, Frequencies* freq) const {
		double E = potential::totalEnergy(*pot, point);
		if (E > 0) {
			printf("Energy is: %f and positive so no orbit\n", E);
			return ActionAngles(actions::Actions(NAN, NAN, NAN), actions::Angles(NAN, NAN, NAN));
		}
		const double tol = 1e-8;
		double phi0 = point.phi, last_diff = 1e6;
		while (fabs(phi0) > M_PI) phi0 += phi0 > M_PI ? -M_PI : M_PI;
		coord::PosMomCyl P0(point.R, point.z, phi0, point.vR, point.vz, point.vphi * point.R);
		Angles trueA;
		ActionAngles aa(AF.actionAngles(point));
		Actions J(aa);
		Torus T1(TG.fitTorus(J));
		Angles aT = T1.GF.toyA(aa);
		std::vector<double> rep;
		int kount = 0;
		while (kount < 10) {
			/*if (kount == 1) {
				T = Torus(TG.fitTorus(J));
			}
			else if (kount > 0) T = TG.fitFrom(J, T);*/
			Torus T(TG.fitTorus(J));
			coord::PosMomCyl P(T.from_toy(aT));
			double phi_diff = P0.phi - P.phi;
			if (fabs(phi_diff) > M_PI) {
				P.phi += phi_diff > 0 ? 2 * M_PI : -2 * M_PI;
			}
			/*std::vector<double> dRp = {P0.R - P.R,P0.z - P.z,P0.phi - P.phi,
			P0.pR - P.pR,P0.pz - P.pz,P0.pphi - P.pphi };
			//*/
			
			std::vector<double> dRp = { P0.R - P.R,P0.z - P.z,
			P0.pR - P.pR,P0.pz - P.pz };
			double diff = 0;
			//for (int i = 0; i < 6; i++) diff += pow_2(dRp[i]);
			for (int i = 0; i < 4; i++) diff += pow_2(dRp[i]);
			rep.push_back(diff);
			//
			if ((diff)< tol) {
				aT.thetaphi = T.TM.pq2aa(P0).thetaphi;
				//aT.thetaphi = T.Is.pq2aa(P0).thetaphi;
				trueA = T.GF.trueA(aT);
				if (freq) *freq = T.freqs;
				break;
			}
			DerivAngCyl dRdt; DerivAct<coord::Cyl> dRdJ;
			ActionAngles AA(J, aT);
			T.PosDerivs(AA, dRdt); T.PosDerivJ(AA, dRdJ);
			//math::Matrix<double> M(6, 6);
			/*
			M(0, 0) = dRdt.dbythetar.R; M(0, 1) = dRdt.dbythetaz.R; M(0, 2) = dRdt.dbythetaphi.R;
			M(0, 3) = dRdJ.dbyJr.R; M(0, 4) = dRdJ.dbyJz.R; M(0, 5) = dRdJ.dbyJphi.R;
			M(1, 0) = dRdt.dbythetar.z; M(1, 1) = dRdt.dbythetaz.z; M(1, 2) = dRdt.dbythetaphi.z;
			M(1, 3) = dRdJ.dbyJr.z; M(1, 4) = dRdJ.dbyJz.z; M(1, 5) = dRdJ.dbyJphi.z;
			M(2, 0) = dRdt.dbythetar.phi; M(2, 1) = dRdt.dbythetaz.phi; M(2, 2) = dRdt.dbythetaphi.phi;
			M(2, 3) = dRdJ.dbyJr.phi; M(2, 4) = dRdJ.dbyJz.phi; M(2, 5) = dRdJ.dbyJphi.phi;
			M(3, 0) = dRdt.dbythetar.pR; M(3, 1) = dRdt.dbythetaz.pR; M(3, 2) = dRdt.dbythetaphi.pR;
			M(3, 3) = dRdJ.dbyJr.pR; M(3, 4) = dRdJ.dbyJz.pR; M(3, 5) = dRdJ.dbyJphi.pR;
			M(4, 0) = dRdt.dbythetar.pz; M(4, 1) = dRdt.dbythetaz.pz; M(4, 2) = dRdt.dbythetaphi.pz;
			M(4, 3) = dRdJ.dbyJr.pz; M(4, 4) = dRdJ.dbyJz.pz; M(4, 5) = dRdJ.dbyJphi.pz;
			M(5, 0) = dRdt.dbythetar.pphi; M(5, 1) = dRdt.dbythetaz.pphi; M(5, 2) = dRdt.dbythetaphi.pphi;
			M(5, 3) = dRdJ.dbyJr.pphi; M(5, 4) = dRdJ.dbyJz.pphi; M(5, 5) = dRdJ.dbyJphi.pphi;
			*/
			math::Matrix<double> M(4, 4);
			M(0, 0) = dRdt.dbythetar.R; M(0, 1) = dRdt.dbythetaz.R; 
			M(0, 2) = dRdJ.dbyJr.R; M(0, 3) = dRdJ.dbyJz.R; 
			M(1, 0) = dRdt.dbythetar.z; M(1, 1) = dRdt.dbythetaz.z; 
			M(1, 2) = dRdJ.dbyJr.z; M(1, 3) = dRdJ.dbyJz.z;
			M(2, 0) = dRdt.dbythetar.pR; M(2, 1) = dRdt.dbythetaz.pR;
			M(2, 2) = dRdJ.dbyJr.pR; M(2, 3) = dRdJ.dbyJz.pR; 
			M(3, 0) = dRdt.dbythetar.pz; M(3, 1) = dRdt.dbythetaz.pz;
			M(3, 2) = dRdJ.dbyJr.pz; M(3, 3) = dRdJ.dbyJz.pz;
			
			

			math::LUDecomp LUM(M);
			std::vector<double> dJt = LUM.solve(dRp);
			aT.thetar += dJt[0]; aT.thetaz += dJt[1];
			//aT.thetaphi += dJt[2];
			//J.Jr += dJt[3]; J.Jz += dJt[4];
			//J.Jphi += dJt[5];
			if (fabs(dJt[2]) > 0.2 * J.Jr) {
							J.Jr *= (1 + 0.2 * math::sign(dJt[2]));
						}
						else {
							J.Jr += dJt[2];
						}if (fabs(dJt[3]) > 0.2 * J.Jz) {
							J.Jz *= (1 + 0.2 * math::sign(dJt[3]));
						}
						else {
							J.Jz += dJt[3];
			}
			//*/

			last_diff = diff; kount++;
		}
		//std::cout << kount <<" "<<diff1<< '\n';
		if (kount == 9) {
			printf("failed to converge\n");
		}
		/*if (kount > 0) {
			//printf("kount %d ", kount);
			for (int i = 0;i < rep.size();i++) printf("%g ", rep[i]);
			printf("\n");
		}*/
		return ActionAngles(J, trueA);
	}
	//end fromT
}//namespace
