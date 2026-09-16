#include <octave/oct.h>

DEFUN_DLD (loopPrimitivesNuc, args, nargout,
           "This function takes as arguments basis_a,basis_b,Z,AL, and calculates the vectors"
           "xValues,KabValues,RPAValues,RPBValues,RPCValues,pValues,WeightValues,ZValues")
{

  Matrix basis_a = args(0).matrix_value ();
  Matrix basis_b = args(1).matrix_value ();
  Matrix Z       = args(2).matrix_value ();
  Matrix AL      = args(3).matrix_value ();


  //octave_stdout << "loopPrimitivesNuc "
  //              << args.length() << " input arguments and "
  //              << nargout << " output arguments.\n";

  //octave_idx_type na;
  //octave_idx_type nb;
  //octave_idx_type N;

  //octave_stdout << basis_a(0) << "\n";
  //octave_stdout << basis_b(0) << "\n";
  //octave_stdout << AL.rows() << "\n";
  const int Natoms = AL.rows();
  const int Nints = basis_a(0)*basis_b(0)*AL.rows();

  //octave_stdout << Nints << "\n";

  dim_vector dv1 (Nints,1);
  dim_vector dv3 (Nints,3);

  Matrix xValues (dv1);
  Matrix KabValues (dv1);
  Matrix RPAValues (dv3);
  Matrix RPBValues (dv3);
  Matrix RPCValues (dv3);
  Matrix pValues (dv1);
  Matrix WeightValues (dv1);
  Matrix ZValues (dv1);

  double xa = basis_a(2);
  double ya = basis_a(3);
  double za = basis_a(4);

  double xb = basis_b(2);
  double yb = basis_b(3);
  double zb = basis_b(4);

  Matrix RAB (3,1);
  RAB(0) = xa-xb;
  RAB(1) = ya-yb;
  RAB(2) = za-zb; //column vector

  int t = 0;

  for (octave_idx_type na = 0; na < basis_a(0); na++){
      int La = basis_a(1);
      int Dima = (La+1)*(La+2)/2;
      double alphaa = basis_a(5+Dima+(na)*3+0);
      double ca = basis_a(5+Dima+(na)*3+1);
      double Na = basis_a(5+Dima+(na)*3+2);


      for (octave_idx_type nb = 0; nb < basis_b(0); nb++){
          int Lb = basis_b(1);
          int Dimb = (Lb+1)*(Lb+2)/2;
          double alphab = basis_b(5+Dimb+(nb)*3+0);
          double cb = basis_b(5+Dimb+(nb)*3+1);
          double Nb = basis_b(5+Dimb+(nb)*3+2);

          double p = alphaa+alphab;

          double Px = (alphaa*xa + alphab*xb)/p;
          double Py = (alphaa*ya + alphab*yb)/p;
          double Pz = (alphaa*za + alphab*zb)/p;

          double PAx = Px-xa;
          double PAy = Py-ya;
          double PAz = Pz-za;

          double PBx = Px-xb;
          double PBy = Py-yb;
          double PBz = Pz-zb;

          double rhoAB = alphaa*alphab/p;
          double Kab = exp(-rhoAB*(pow(RAB(0),2)+pow(RAB(1),2)+pow(RAB(2),2)));
          //octave_stdout << alphaa << "\n";
          //octave_stdout << alphab << "\n";
          //octave_stdout << rhoAB << "\n";
          //octave_stdout << Kab << "\n";


          for (octave_idx_type iatom = 0; iatom < Natoms; iatom++){

              double PCx = Px-AL(iatom,0);
              double PCy = Py-AL(iatom,1);
              double PCz = Pz-AL(iatom,2);
              double RPC2 = pow(PCx,2)+pow(PCy,2)+pow(PCz,2);

              double x = p*RPC2;

              ZValues(t) = Z(iatom);
              xValues(t) = x;
              KabValues(t) = Kab;
              RPAValues(t,0) = PAx;
              RPAValues(t,1) = PAy;
              RPAValues(t,2) = PAz;
              RPBValues(t,0) = PBx;
              RPBValues(t,1) = PBy;
              RPBValues(t,2) = PBz;
              RPCValues(t,0) = PCx;
              RPCValues(t,1) = PCy;
              RPCValues(t,2) = PCz;
              pValues(t) = p;
              WeightValues(t) = ca*Na*cb*Nb;

              t = t + 1;


          }

      }
      //octave_stdout << na << "\n";
  }


  octave_value_list retval (nargout);

  retval(0) = octave_value(xValues);
  retval(1) = octave_value(KabValues);
  retval(2) = octave_value(RPAValues);
  retval(3) = octave_value(RPBValues);
  retval(4) = octave_value(RPCValues);
  retval(5) = octave_value(pValues);
  retval(6) = octave_value(WeightValues);
  retval(7) = octave_value(ZValues);
  return retval;
}
