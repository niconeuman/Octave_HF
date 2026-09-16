#include <octave/oct.h>

DEFUN_DLD (loopPrimitives, args, nargout,
           "This function takes as arguments pair_data2_ab,pair_data2_cd, and calculates the vectors"
           "RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,ooppqValues,xValues,PrefactorValues,WeightValues")
{

    Matrix pair_data2_ab = args(0).matrix_value ();
    Matrix pair_data2_cd = args(1).matrix_value ();

    const int Nintsab = pair_data2_ab.rows();
    const int Nintscd = pair_data2_cd.rows();
    const int Nints = pair_data2_ab.rows()*pair_data2_cd.rows();

    //octave_stdout << Nints << "\n";

    dim_vector dv1 (Nints,1);
    dim_vector dv3 (Nints,3);

    Matrix xValues (dv1);
    Matrix KabValues (dv1);
    Matrix RPAValues (dv3);
    Matrix RPBValues (dv3);
    Matrix RQCValues (dv3);
    Matrix RQDValues (dv3);
    Matrix RWPValues (dv3);
    Matrix RWQValues (dv3);
    Matrix pValues (dv1);
    Matrix qValues (dv1);
    Matrix ppqValues (dv1);
    Matrix ooppqValues (dv1);
    Matrix WeightValues (dv1);
    Matrix PrefactorValues (dv1);

    int t = 0;

    for (octave_idx_type nab = 0; nab < Nintsab; nab++){

        // Matrix pValuesab  = pair_data2_ab(:,0);
        // Matrix PxValuesab = pair_data2_ab(:,1);
        // Matrix PyValuesab = pair_data2_ab(:,2);
        // Matrix PzValuesab = pair_data2_ab(:,3);
        // Matrix KabValues = pair_data2_ab(:,4);
        // Matrix WeightValuesab = pair_data2_ab(:,5);
        // Matrix RPAValuesab = pair_data2_ab(:,6:8);
        // Matrix RPBValuesab = pair_data2_ab(:,9:11);
        // Matrix pPxValuesab = pair_data2_ab(:,12);
        // Matrix pPyValuesab = pair_data2_ab(:,13);
        // Matrix pPzValuesab = pair_data2_ab(:,14);

        double p = pair_data2_ab(nab,0);

        double Px = pair_data2_ab(nab,1);
        double Py = pair_data2_ab(nab,2);
        double Pz = pair_data2_ab(nab,3);

        double Kab = pair_data2_ab(nab,4);
        double weightab = pair_data2_ab(nab,5);

        double PAx = pair_data2_ab(nab,6);
        double PAy = pair_data2_ab(nab,7);
        double PAz = pair_data2_ab(nab,8);

        double PBx = pair_data2_ab(nab,9);
        double PBy = pair_data2_ab(nab,10);
        double PBz = pair_data2_ab(nab,11);

        double pPx = pair_data2_ab(nab,12);
        double pPy = pair_data2_ab(nab,13);
        double pPz = pair_data2_ab(nab,14);


        for (octave_idx_type ncd = 0; ncd < Nintscd; ncd++){

            double q = pair_data2_cd(ncd,0);

            double Qx = pair_data2_cd(ncd,1);
            double Qy = pair_data2_cd(ncd,2);
            double Qz = pair_data2_cd(ncd,3);

            double Kcd = pair_data2_cd(ncd,4);
            double weightcd = pair_data2_cd(ncd,5);

            double QCx = pair_data2_cd(ncd,6);
            double QCy = pair_data2_cd(ncd,7);
            double QCz = pair_data2_cd(ncd,8);

            double QDx = pair_data2_cd(ncd,9);
            double QDy = pair_data2_cd(ncd,10);
            double QDz = pair_data2_cd(ncd,11);

            double qQx = pair_data2_cd(ncd,12);
            double qQy = pair_data2_cd(ncd,13);
            double qQz = pair_data2_cd(ncd,14);

            //Now I start calculating quartet quantities
            double PQx = Px-Qx;
            double PQy = Py-Qy;
            double PQz = Pz-Qz;

            double PQ2 = PQx*PQx+PQy*PQy+PQz*PQz;

            double ppq = p+q;
            double ooppq = 1/ppq;
            double alpha = (p*q)*ooppq;
            double x = alpha*PQ2;

            double Wx = (pPx+qQx)*ooppq;
            double Wy = (pPy+qQy)*ooppq;
            double Wz = (pPz+qQz)*ooppq;

            double WPx = Wx-Px;
            double WPy = Wy-Py;
            double WPz = Wz-Pz;

            double WQx = Wx-Qx;
            double WQy = Wy-Qy;
            double WQz = Wz-Qz;

            double Prefactor = Kab*Kcd*2*17.493418327624862/(p*q*sqrt(ppq));

            xValues(t) = x;

            RPAValues(t,0) = PAx;
            RPAValues(t,1) = PAy;
            RPAValues(t,2) = PAz;

            RPBValues(t,0) = PBx;
            RPBValues(t,1) = PBy;
            RPBValues(t,2) = PBz;

            RQCValues(t,0) = QCx;
            RQCValues(t,1) = QCy;
            RQCValues(t,2) = QCz;

            RQDValues(t,0) = QDx;
            RQDValues(t,1) = QDy;
            RQDValues(t,2) = QDz;

            pValues(t) = p;
            qValues(t) = q;
            ppqValues(t) = ppq;
            ooppqValues(t) = ooppq;

            RWPValues(t,0) = WPx;
            RWPValues(t,1) = WPy;
            RWPValues(t,2) = WPz;

            RWQValues(t,0) = WQx;
            RWQValues(t,1) = WQy;
            RWQValues(t,2) = WQz;

            WeightValues(t) = weightab*weightcd;
            PrefactorValues(t) = Prefactor;
            t = t + 1;

        }
        //octave_stdout << na << "\n";
    }
    //RPAValues,RPBValues,RQCValues,RQDValues,RWPValues,RWQValues,pValues,qValues,ppqValues,ooppqValues,xValues,PrefactorValues,WeightValues
    //
    octave_value_list retval (nargout);

    retval(0) = octave_value(RPAValues);
    retval(1) = octave_value(RPBValues);

    retval(2) = octave_value(RQCValues);
    retval(3) = octave_value(RQDValues);

    retval(4) = octave_value(RWPValues);
    retval(5) = octave_value(RWQValues);

    retval(6) = octave_value(pValues);
    retval(7) = octave_value(qValues);
    retval(8) = octave_value(ppqValues);
    retval(9) = octave_value(ooppqValues);

    retval(10) = octave_value(xValues);

    retval(11) = octave_value(PrefactorValues);
    retval(12) = octave_value(WeightValues);

    return retval;
}
