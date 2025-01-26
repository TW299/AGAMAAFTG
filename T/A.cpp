#include "actions_staeckel.h"
#include "orbit.h"
#include "actions_newtorus.h"
#include "potential_factory.h"
#include "units.h"

#include <vector>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <fstream>
#include <stdio.h>
const units::InternalUnits unit(units::galactic_Myr);
/*To get torus code you first need to get a pointer to the potential with ptrpotential. Then you need to initialise a Torus Generator from that pterpotential 
and action finder. Then intiialise the class actionfinderTG, and AFTG.actionagnles(xv) where xv is posvelcyl coordinates wil get angle action coordinates.
The code:

actions::Frequencies f;
actions::ActionAngles aa=AFTG.actionAngles(xv,&f);

will get the frequencies as a class f (storing f.omegar,f.omegaz,f.omegaphi) and action angles aa.
  */
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
 //Torus map J and theta to some momentum position.
    coord::PosMomCyl xp = T.from_true(theta);
    
    coord::PosVelCyl xv(xp.R, xp.z, xp.phi, xp.pR, xp.pz, xp.pphi / xp.R);

    const actions::ActionFinderAxisymFudge AF(pot);
    actions::Frequencies f = T.freqs;
    int N = 200;
    double duration = 2 * M_PI * 2 / f.Omegar;
    std::cout << duration << '\n';
    double dt = duration/N;
    std::vector<std::pair<coord::PosVelCyl, double> > traj(orbit::integrateTraj(xv,duration,dt,*pot));
    std::vector<double> JTGr(traj.size(), 0.0);std::vector<double> Jrs(traj.size(), 0.0);
    std::ofstream file;
    actions::ActionFinderTG AFTG(pot, AF, TG);
  //rewrites or writes new text file in the diretory titled File6.dat containing angle and action information.
    file.open("File6.dat");
    file << "Jr Jz thetar thetaz Jrs Jzs time"<<'\n';
  
    for (int i = 0;i < traj.size();i++) {
        coord::PosVelCyl pv1 = traj[i].first;
        
        actions::ActionAngles aa1 = AFTG.actionAngles(pv1);
        actions::Actions JS = AF.actions(pv1);
        file << aa1.Jr << " " << aa1.Jz << " "<<aa1.thetar<<" "<<aa1.thetaz<<" " << JS.Jr << " " << JS.Jz  << " " <<traj[i].second << '\n';
    }
    file.close();
    return 0;
}
