#include "actions_newisochrone.h"

namespace actions {

	namespace {

		class Is_Horse {
		private:
			double Js, b;
			double Js2, Js4, b2, c, H, Omegar, Omegaphi, Jr, Jz, Jphi, thetar, thetaz, thetaphi, L, Ls, frat;
			double a, ap, boc, u, e, sini, cosi, eta, cseta, sneta, snheta, csheta, tnheta, psi, snpsi, cspsi;
			double r, pr, pth, csth, snth;  int signJphi;
			void derivs1(double* dlncdJ, double* dedJ, double* detadJ,
				double* dudJ, double* dpsidJ, double* dchidJ);
			void derivs2(double* dlncdJ, double* dedJ, double* detadJ, double* dudJ,
				double* dpsidJ, double* dchidJ, double* drdJ, double* dthdJ,
				double* dphidJ, double* dprdJ, double* dpthdJ);
		public:
			Is_Horse(const double& _Js, const double& _b) :
				Js(_Js), Js2(_Js* _Js), Js4(Js2* Js2), b(_b), b2(_b* _b) {
			}
			coord::PosMomSph aa2pq(const ActionAngles&, Frequencies*,
				DerivAct<coord::Sph>*, DerivAng<coord::Sph>*);
			ActionAngles pq2aa(const coord::PosMomSph&, Frequencies* freqs);
		};
		void Is_Horse::derivs1(double* dlncdJ, double* dedJ, double* detadJ,
			double* dudJ, double* dpsidJ, double* dchidJ) {
			double X = .5 * Js2 / (pow_2(H) * b * c);
			dlncdJ[0] = X * Omegar; dlncdJ[1] = X * Omegaphi; dlncdJ[2] = X * Omegaphi;
			X = .5 * (1 / e - e);
			dedJ[0] = X * (1 + boc / (1 + boc)) * dlncdJ[0];
			dedJ[1] = X * (1 + boc / (1 + boc)) * dlncdJ[1] - 2 * X / L; dedJ[2] = dedJ[1];
			//X=(eta-thetar)*(1+boc)/(u+boc);
			X = sneta / (u + boc);//avoid division by e
			for (int i = 0; i < 3; i++) {
				detadJ[i] = X * (dedJ[i] + e * boc / (1 + boc) * dlncdJ[i]);
				dudJ[i] = e * sneta * detadJ[i] - cseta * dedJ[i];
			}
			X = a / (1 - e * e);
			double Y = ap / (pow_2(1 + 2 * boc) - e * e);
			double dadJ[3], dapdJ[3];
			for (int i = 0; i < 3; i++) {
				dadJ[i] = X * dedJ[i];
				dapdJ[i] = Y * ((1 + 2 * boc) * dedJ[i] - 2 * e * boc * dlncdJ[i]);
			}
			X = 1 / (1 + pow_2(a * tnheta)); Y = .5 * a / pow_2(csheta);
			double dfadJ[3], dfapdJ[3];
			for (int i = 0; i < 3; i++)
				dfadJ[i] = X * (tnheta * dadJ[i] + Y * detadJ[i]);
			X = 1 / (1 + pow_2(ap * tnheta)); Y = .5 * ap / pow_2(csheta);
			for (int i = 0; i < 3; i++)
				dfapdJ[i] = X * (tnheta * dapdJ[i] + Y * detadJ[i]);
			double dfratdL = 2 * pow_2(Js / L) / pow(1 + 4 * pow_2(Js / L), 1.5) / L;
			X = 2 * atan(ap * tnheta) - thetar; Y = 2 * frat - 1;
			for (int i = 0; i < 3; i++)
				dpsidJ[i] = dfadJ[i] + Y * dfapdJ[i];
			for (int i = 1;i < 3;i++)
				dpsidJ[i] += X * dfratdL;
			dchidJ[0] = fabs(Jphi) / L * dpsidJ[0];
			dchidJ[1] = -fabs(Jphi) / pow_2(L) * snpsi * cspsi + fabs(Jphi) / L * dpsidJ[1];
			dchidJ[2] = (signJphi / L - fabs(Jphi) / pow_2(L)) * snpsi * cspsi + fabs(Jphi) / L * dpsidJ[2];
		}
		void Is_Horse::derivs2(double* dlncdJ, double* dedJ, double* detadJ,
			double* dudJ, double* dpsidJ, double* dchidJ,
			double* drdJ, double* dthdJ, double* dphidJ,
			double* dprdJ, double* dpthdJ) {
			double tnpsi = snpsi / cspsi;
			dthdJ[0] = -sini * cspsi / snth * dpsidJ[0];
			dthdJ[1] = -snpsi / (L * sini * snth) * pow_2(Jphi / L)
				- sini * cspsi / snth * dpsidJ[1];
			dthdJ[2] = snpsi / (L * sini * snth) * (1 - Jphi / L) * (Jphi / L)
				- sini * cspsi / snth * dpsidJ[2];
			dpthdJ[0] = L * sini * (snpsi * dpsidJ[0] / snth + cspsi * csth * dthdJ[0] / pow_2(snth));
			dpthdJ[1] = -cspsi / (snth * sini)
				+ L * sini * (snpsi * dpsidJ[1] / snth + cspsi * csth * dthdJ[1] / pow_2(snth));
			if (std::isnan(dpthdJ[0])) {
				printf("TM error %f %f %f %f %f %f %f %f\n", pth, L, Jphi, tnpsi, dpsidJ[1], csth, snth, dthdJ[1]);
				exit(0);
			}
			dpthdJ[2] = -(1 - Jphi / L) * cspsi / (snth * sini)
				+ L * sini * (snpsi * dpsidJ[2] / snth + cspsi * csth * dthdJ[2] / pow_2(snth));
			for (int i = 0; i < 3; i++) {
				drdJ[i] = c * c / r * ((u + boc) * dudJ[i] - u * boc * dlncdJ[i]) + r * dlncdJ[i];
				dprdJ[i] = .5 * pr * (dlncdJ[i] + boc / (1 + boc) * dlncdJ[i] + 2 / e * dedJ[i] - 2 / r * drdJ[i])
					+ Js * sqrt(1 / (boc * (1 + boc))) * e / r * cseta * detadJ[i];
				dphidJ[i] = signJphi * dchidJ[i];
			}
		}
		coord::PosMomSph Is_Horse::aa2pq(const ActionAngles& aa, Frequencies* freqs,
			DerivAct<coord::Sph>* dJ, DerivAng<coord::Sph>* dA) {
			Jr = aa.Jr; Jz = aa.Jz; Jphi = aa.Jphi;
			thetar = math::wrapAngle(aa.thetar);
			if (thetar > M_PI) thetar -= 2 * M_PI;
			thetaz = math::wrapAngle(aa.thetaz);
			thetaphi = math::wrapAngle(aa.thetaphi);
			signJphi = Jphi > 0 ? 1 : -1;
			L = Jz + fabs(Jphi);
			cosi = Jphi / L; sini = sqrt(1 - cosi * cosi);
			Ls = sqrt(L * L + 4 * Js2);
			double bot = aa.Jr + .5 * (L + Ls);
			H = -.5 * Js4 / b2 / pow_2(bot);
			Omegar = Js4 / b2 / pow(bot, 3);
			frat = .5 * (1 + L / Ls);
			Omegaphi = frat * Omegar;
			if (freqs) {
				freqs->Omegar = Omegar;
				freqs->Omegaz = Omegaphi;
				freqs->Omegaphi = Omegaphi;
			}
			c = Js2 / (-2 * b * H) - b;
			boc = b / c;
			const double e2 = 1 - L * L * boc / Js2 * (1 + boc);
			e = e2 > 0 ? sqrt(e2) : 1e-20;// e=0 is problematic
			eta = math::solveKepler(e / (1 + boc), thetar);
			math::sincos(eta, sneta, cseta);
			math::sincos(.5 * eta, snheta, csheta); tnheta = snheta / csheta;
			u = 1 - e * cseta;
			a = sqrt((1 + e) / (1 - e));
			ap = sqrt((1 + e + 2 * boc) / (1 - e + 2 * boc));
			//the following eq is only valid for -pi<=thetar<pi
			psi = thetaz - frat * thetar + atan(a * tnheta) + (2 * frat - 1) * atan(ap * tnheta);
			math::sincos(psi, snpsi, cspsi);
			r = c * sqrt(u * (u + 2 * boc));
			csth = sini * snpsi;
			snth = sqrt(1 - csth * csth);
			double chi = Jz != 0 ? atan2(fabs(Jphi) * snpsi, L * cspsi) : psi;
			double phi = math::wrapAngle(thetaphi + (chi - thetaz) * signJphi);
			pr = Js / sqrt(boc * (1 + boc)) * e / r * sneta;
			pth = -L * sini * cspsi / snth;
			coord::PosMomSph p(r, acos(csth), phi, pr, pth, aa.Jphi);
			if (dJ) {// Compute d/dJi
				double dlncdJ[3], dedJ[3], detadJ[3], dudJ[3], dpsidJ[3], dchidJ[3];
				derivs1(dlncdJ, dedJ, detadJ, dudJ, dpsidJ, dchidJ);
				double drdJ[3], dthdJ[3], dphidJ[3], dprdJ[3], dpthdJ[3];
				derivs2(dlncdJ, dedJ, detadJ, dudJ, dpsidJ, dchidJ,
					drdJ, dthdJ, dphidJ, dprdJ, dpthdJ);
				dJ->dbyJr.r = drdJ[0];    dJ->dbyJr.theta = dthdJ[0];
				dJ->dbyJz.r = drdJ[1];    dJ->dbyJz.theta = dthdJ[1];
				dJ->dbyJphi.r = drdJ[2]; dJ->dbyJphi.theta = dthdJ[2];
				dJ->dbyJr.phi = dphidJ[0]; dJ->dbyJz.phi = dphidJ[1]; dJ->dbyJphi.phi = dphidJ[2];
				dJ->dbyJr.pr = dprdJ[0]; dJ->dbyJr.ptheta = dpthdJ[0];
				dJ->dbyJz.pr = dprdJ[1]; dJ->dbyJz.ptheta = dpthdJ[1];
				dJ->dbyJphi.pr = dprdJ[2];  dJ->dbyJphi.ptheta = dpthdJ[2];
			} if (dA) {// compute d/dtheta
				double dfadeta = .5 * a / (pow_2(csheta) + pow_2(a * snheta));
				double dfapdeta = .5 * ap / (pow_2(csheta) + pow_2(ap * snheta));
				double detadthr = (1 + boc) / (u + boc);
				double dpsidthr = -frat + (dfadeta + (2 * frat - 1) * dfapdeta) * detadthr;
				double dpsidthz = 1;
				double cmn = fabs(Jphi) / (L * (pow_2(cspsi) + pow_2(Jphi / L * snpsi)));
				double dchidthr = cmn * dpsidthr;
				double dchidthz = cmn * dpsidthz;

				dA->dbythetar.r = c * c / r * (1 + boc) * e * sneta;
				dA->dbythetaz.r = 0;
				dA->dbythetaphi.r = 0;
				dA->dbythetar.theta = pth / L * dpsidthr;
				dA->dbythetaz.theta = pth / L;
				dA->dbythetaphi.theta = 0;
				dA->dbythetar.phi = signJphi * dchidthr;
				dA->dbythetaz.phi = signJphi * (dchidthz - 1);
				dA->dbythetaphi.phi = 1;
				dA->dbythetar.pr = -pr / r * dA->dbythetar.r
					+ Js / sqrt(boc * (1 + boc)) * e / r * cseta * detadthr;
				dA->dbythetaz.pr = 0;
				dA->dbythetaz.ptheta = (L - pow_2(pth) / L) * csth / snth;
				dA->dbythetar.ptheta = dA->dbythetaz.ptheta * dpsidthr;
				dA->dbythetaphi.pr = dA->dbythetaphi.ptheta = 0;
				dA->dbythetaphi.pphi = 0;
			}
			return p;
		}
		ActionAngles Is_Horse::pq2aa(const coord::PosMomSph& p, Frequencies* freqs) {
			math::sincos(p.theta, snth, csth);
			L = p.pphi == 0 ? fabs(p.ptheta) : sqrt(pow_2(p.ptheta) + pow_2(p.pphi / snth));
			Ls = sqrt(L * L + 4 * Js2);
			Jz = L - fabs(p.pphi);
			H = .5 * (pow_2(p.pr) + pow_2(p.ptheta / p.r)) - Js2 / b / (b + sqrt(b * b + p.r * p.r));
			if (p.pphi != 0) H += .5 * pow_2(p.pphi / (p.r * snth));
			double rtH = sqrt(-2 * H);
			Jr = Js2 / (b * rtH) - .5 * (L + Ls);
			double bot = Jr + .5 * (L + Ls);
			Omegar = Js4 / b2 / pow(bot, 3);
			frat = .5 * (1 + L / Ls);
			Omegaphi = frat * Omegar;
			if (freqs) {
				freqs->Omegar = Omegar;
				freqs->Omegaz = Omegaphi;
				freqs->Omegaphi = Omegaphi;
			}
			c = Js2 / (-2 * b * H) - b;
			boc = b / c;
			const double e2 = 1 - L * L * boc / Js2 * (1 + boc);
			e = e2 > 0 ? sqrt(e2) : 1e-20;
			eta = atan2(p.r * p.pr / rtH, b + c - sqrt(b * b + p.r * p.r));
			psi = atan2(csth, -p.ptheta / L * snth);
			thetar = eta - e / (1 + boc) * sin(eta);
			tnheta = tan(.5 * eta);
			a = sqrt((1 + e) / (1 - e)); ap = sqrt((1 + e + 2 * boc) / (1 - e + 2 * boc));
			thetaz = math::wrapAngle(psi + frat * thetar - atan(a * tnheta) - (2 * frat - 1) * atan(ap * tnheta));
			double chi = atan2(csth / snth * p.pphi, -p.ptheta);
			thetaphi = math::wrapAngle(p.phi - chi + math::sign(p.pphi) * thetaz);
			return 	ActionAngles(Actions(Jr, Jz, p.pphi), Angles(thetar, thetaz, thetaphi));
		}

	}//namespace anon

	double Isochrone::H(coord::PosMomSph& rp) const {
		double H = .5 * (pow_2(rp.pr) + pow_2(rp.ptheta / rp.r))
			- Js * Js / b / (b + sqrt(b * b + pow_2(rp.r)));
		return rp.pphi == 0 ? H : H + .5 * pow_2(rp.pphi / (rp.r * sin(rp.theta)));
	}
	double Isochrone::PE(coord::PosMomSph& rp) const {
		return -Js * Js / b / (b + sqrt(b * b + pow_2(rp.r)));
	}
	Actions Isochrone::pq2J(const coord::PosMomSph& p) const {
		double snt, cst, Js2 = Js * Js;
		math::sincos(p.theta, snt, cst);
		double L = p.pphi == 0 ? fabs(p.ptheta) : sqrt(pow_2(p.ptheta) + pow_2(p.pphi / snt));
		double Ls = sqrt(L * L + 4 * Js2);
		double Jz = L - fabs(p.pphi);
		double H = .5 * (pow_2(p.pr) + pow_2(p.ptheta / p.r)) - Js2 / b / (b + sqrt(b * b + p.r * p.r));
		if (p.pphi != 0) H += .5 * pow_2(p.pphi / (p.r * snt));
		double rtH = sqrt(-2 * H);
		double Jr = Js2 / (b * rtH) - .5 * (L + Ls);
		return Actions(Jr, Jz, p.pphi);
	}
	coord::PosMomSph Isochrone::aa2pq(const ActionAngles& aa, Frequencies* freqs,
		DerivAct<coord::Sph>* dJ, DerivAng<coord::Sph>* dA) const {
		Is_Horse WH(Js, b);
		return WH.aa2pq(aa, freqs, dJ, dA);
	}
	ActionAngles Isochrone::pq2aa(const coord::PosMomSph& rtheta, Frequencies* freqs) const {
		Is_Horse WH(Js, b);
		return WH.pq2aa(rtheta, freqs);
	}

}//namespace actions