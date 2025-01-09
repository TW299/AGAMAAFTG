#include "actions_staeckel.h"
#include "orbit.h"
#include <vector>
#include <iomanip>
#include "actions_newtorus.h"
#include "potential_composite.h"
#include "potential_factory.h"
#include "potential_analytic.h"
#include "units.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <stdio.h>const units::InternalUnits unit(units::galactic_Myr);

int main() {
  //NFW potential with q=0.5
    potential::PtrPotential pot = potential::createPotential(utils::KeyValueMap("type=spheroid, gamma=1, beta=3,alpha=1, scaleradius=1,densityNorm=0.1, q=0.5"));
  //second variable is tolerance of J
    const actions::TorusGenerator TG(*pot, 1e-6);
    double Jr = 0.1, Jz = 0.55, Jphi = 3.0;
    actions::Actions J(Jr, Jz, Jphi);
    double ar = 2, az = .5, arp = -ar, azp = M_PI - az;
    actions::Torus T(TG.fitTorus(J));
    actions::Angles theta(1.0, 1.3, 1.3);
    coord::PosMomCyl xp = T.from_true(theta);
    
    coord::PosVelCyl xv(xp.R, xp.z, xp.phi, xp.pR, xp.pz, xp.pphi / xp.R);

    const actions::ActionFinderAxisymFudge AF(pot);
    actions::Frequencies f = T.freqs;
    int N = 200;
    double duration = 2 * M_PI * 2 / f.Omegar;
    std::cout << duration << '\n';
    double dt = duration/N;
    std::vector<std::pair<coord::PosVelCyl, double> > traj(T.orbit(theta, dt, duration));
    std::vector<double> JTGr(traj.size(), 0.0);std::vector<double> Jrs(traj.size(), 0.0);
    std::ofstream file;
    actions::ActionFinderTG AFTG(pot, AF, TG);
    file.open("File6.dat");
    file << "Jr Jz thetar thetaz Jrs Jzs time"<<'\n';
    for (int i = 0;i < traj.size();i++) {
        coord::PosVelCyl pv1 = traj[i].first;
        
        actions::ActionAngles aa1 = AFTG.actionAngles(pv1);
        actions::Actions JS = AF.actions(pv1);
        std::cout << i << '\n';
        file << aa1.Jr << " " << aa1.Jz << " "<<aa1.thetar<<" "<<aa1.thetaz<<" " << JS.Jr << " " << JS.Jz  << " " <<traj[i].second << '\n';
    }
    file.close();
    return 0;
}
