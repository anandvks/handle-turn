/*
 * Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * Author: Ugo Pattacini
 * email:  ugo.pattacini@iit.it
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * A copy of the license can be found at
 * http://www.robotcub.org/icub/license/gpl.txt
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
*/

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <yarp/os/all.h>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>
#include <cer_kinematics/arm.h>

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace cer::kinematics;


/****************************************************************/
string formatOutput(ArmSolver &s, const int i, const Matrix &Hd,
                    const Vector &q)
{
    Matrix H;
    s.fkin(q,H);
    
    Vector ud=dcm2axis(Hd);
    Vector u=dcm2axis(H);
    ud*=ud[3]; ud.pop_back();
    u*=u[3]; u.pop_back();

    double e_x=norm(Hd.getCol(3).subVector(0,2)-H.getCol(3).subVector(0,2));
    double e_u=norm(ud-u);

    Property p;
    p.put("point",i);

    Bottle b;
    b.addList().read(const_cast<Vector&>(q));
    p.put("q",b.get(0));

    p.put("e_x",e_x);
    p.put("e_u",e_u);

    return p.toString();
}


/****************************************************************/
int main(int argc, char *argv[])
{
    ResourceFinder rf;
    rf.configure(argc,argv);

    if (rf.check("help"))
    {
        cout<<"Options:"<<endl;
        cout<<"--type left|right"<<endl;
        cout<<"--R1-initial-joints \"(0.0 1.0 ... 11.0)\""<<endl;
        cout<<"--handle-length <double> [m]"<<endl;
        cout<<"--rotation <double> [deg]"<<endl;
        cout<<"--points <int>"<<endl;
        cout<<"--output-file <file-name>"<<endl;
        return 0;
    }

    string type=rf.check("type",Value("left")).asString();
    transform(type.begin(),type.end(),type.begin(),::tolower);

    if ((type!="left") && (type!="right"))
    {
        cerr<<"unrecognized type \""<<type<<"\""<<endl;
        return 1;
    }

    double handle_length=rf.check("handle-length",Value(0.08)).asDouble();
    double rotation=rf.check("rotation",Value(10.0)).asDouble();
    int points=rf.check("points",Value(1)).asInt();
    string output_file=rf.check("output-file",Value("output.txt")).asString();

    Vector q0(12);
    //default torso tripod
    q0[0]=q0[1]=q0[2]=0.1343;
    // default torso yaw
    q0[3]=0.0;
    // default arm
    q0[4]=50.509; q0[5]=23.225; q0[6]=-10.107; q0[7]=54.849; q0[8]=-62.853;
    // default wrist tripod
    q0[9]=0.022686; q0[10]=0.020822; q0[11]=0.024994;
    if (Bottle *b=rf.find("R1-initial-joints").asList())
    {
        size_t len=std::min(q0.length(),(size_t)b->size());
        for (size_t i=0; i<len; i++)
            q0[i]=b->get(i).asDouble();
    }

    ArmParameters armp(type);
    ArmSolver solver(armp);

    SolverParameters p=solver.getSolverParameters();
    p.setMode("full_pose+no_torso_heave");
    p.max_cpu_time=1.0;
    solver.setSolverParameters(p);

    Matrix H0;
    solver.fkin(q0,H0);
    
    Matrix C0=eye(4,4);
    C0.setCol(3,H0.getCol(3));
    C0(1,3)+=handle_length;
    
    ofstream fout(output_file.c_str());
    if (!fout.is_open())
    {
        cerr<<"unable to write to \""<<output_file<<"\""<<endl;
        return 2;
    }

    fout<<formatOutput(solver,0,H0,q0)<<endl;

    H0=SE3inv(C0)*H0;
    Vector axis=yarp::math::sign(handle_length)*C0.getCol(0);    

    Vector q=q0;
    for (int i=1; i<=points; i++)
    {
        axis[3]=(M_PI/180.0)*(i*(rotation/points));
        Matrix Hd=C0*axis2dcm(axis)*H0;

        solver.setInitialGuess(q);
        solver.ikin(Hd,q);

        fout<<formatOutput(solver,i,Hd,q)<<endl;
    }

    cout<<"data written to \""<<output_file<<"\""<<endl;
    fout.close();
    return 0;
}


