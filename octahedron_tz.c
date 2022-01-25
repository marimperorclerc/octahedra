#include <math.h>
#include <stdio.h> 
static double
form_volume(double length_a, double b2a_ratio, double c2a_ratio, double t)
{
//octehedron volume formula
// length_a is the half height along the a axis of the octahedron
    return (4./3.) * length_a * (length_a*b2a_ratio) * (length_a*c2a_ratio)*(1.-(1.-t)*(1.-t)*(1.-t));
}

static double
Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double t)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;

   //Integration limits to use in Gaussian quadrature
    const double v1a = 0.0;
    const double v1b = M_PI_2;  //theta integration limits
    const double v2a = 0.0;
    const double v2b = M_PI_2;  //phi integration limits

    double outer_sum = 0.0;
    for(int i=0; i<GAUSS_N; i++) {
        const double theta = 0.5 * ( GAUSS_Z[i]*(v1b-v1a) + v1a + v1b );
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);

        double inner_sum = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            double phi = 0.5 * ( GAUSS_Z[j]*(v2b-v2a) + v2a + v2b );
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            //HERE: Octahedron formula
            const double Qx = q * sin_theta * cos_phi;
    	    const double Qy = q * sin_theta * sin_phi;
    	    const double Qz = q * cos_theta;
    	    const double qx = Qx * length_a;
    	    const double qy = Qy * length_b;
    	    const double qz = Qz * length_c;

            const double pf = 6./(qx*qx-qy*qy)*(1./(1.-(1.-t)*(1.-t)*(1.-t)));
            const double A = -qx*sin(qx)+0.5*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t));
            const double B =  qy*sin(qy)-0.5*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));
            const double AA = A/((qx*qx)-(qz*qz));
            const double BB = B/((qy*qy)-(qz*qz));

            // normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
            const double AP = pf*(AA+BB);


            inner_sum += GAUSS_W[j] * AP * AP;
        }
        inner_sum = 0.5 * (v2b-v2a) * inner_sum;
        outer_sum += GAUSS_W[i] * inner_sum * sin_theta;
    }

    double answer = 0.5*(v1b-v1a)*outer_sum;

    // Normalize by Pi (Eqn. 16).
    // The factor 2 appears because the theta integral has been defined between
    // 0 and pi/2, instead of 0 to pi.
    answer /= M_PI_2; //Form factor P(q)

    // Multiply by contrast^2 and volume^2
    // volume of octahedron
    const double volume = (4./3.)*length_a * length_b * length_c*(1.-(1.-t)*(1.-t)*(1.-t));
    answer *= square((sld-solvent_sld)*volume);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double t)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;

   //Integration limits to use in Gaussian quadrature
    const double v1a = 0.0;
    const double v1b = M_PI_2;  //theta integration limits
    const double v2a = 0.0;
    const double v2b = M_PI_2;  //phi integration limits

    double outer_sum_F1 = 0.0;
    double outer_sum_F2 = 0.0;

    for(int i=0; i<GAUSS_N; i++) {
        const double theta = 0.5 * ( GAUSS_Z[i]*(v1b-v1a) + v1a + v1b );
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);

        double inner_sum_F1 = 0.0;
        double inner_sum_F2 = 0.0;
        for(int j=0; j<GAUSS_N; j++) {
            double phi = 0.5 * ( GAUSS_Z[j]*(v2b-v2a) + v2a + v2b );
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            //HERE: Octahedron formula
            const double Qx = q * sin_theta * cos_phi;
    	    const double Qy = q * sin_theta * sin_phi;
    	    const double Qz = q * cos_theta;
    	    const double qx = Qx * length_a;
    	    const double qy = Qy * length_b;
    	    const double qz = Qz * length_c;
            const double pf = 6./(qx*qx-qy*qy)*(1./(1.-(1.-t)*(1.-t)*(1.-t)));
            const double A = -qx*sin(qx)+0.5*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t));
            const double B =  qy*sin(qy)-0.5*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));
            const double AA = A/((qx*qx)-(qz*qz));
            const double BB = B/((qy*qy)-(qz*qz));

            // normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
            const double AP = pf*(AA+BB);

            inner_sum_F1 += GAUSS_W[j] * AP;
            inner_sum_F2 += GAUSS_W[j] * AP * AP;
        }
        inner_sum_F1 = 0.5 * (v2b-v2a) * inner_sum_F1;
        inner_sum_F2 = 0.5 * (v2b-v2a) * inner_sum_F2;
        outer_sum_F1 += GAUSS_W[i] * inner_sum_F1 * sin_theta;
        outer_sum_F2 += GAUSS_W[i] * inner_sum_F2 * sin_theta;
    }

    outer_sum_F1 *= 0.5*(v1b-v1a);
    outer_sum_F2 *= 0.5*(v1b-v1a);

    // Normalize by Pi (Eqn. 16).
    // The factor 2 appears because the theta integral has been defined between
    // 0 and pi/2, instead of 0 to pi.
    outer_sum_F1 /= M_PI_2;
    outer_sum_F2 /= M_PI_2;

    // Multiply by contrast and volume
    // volume of octahedron
    const double s = (sld-solvent_sld) * (4./3.) *(1.-(1.-t)*(1.-t)*(1.-t))* (length_a * length_b * length_c);

    // Convert from [1e-12 A-1] to [cm-1]
    *F1 = 1e-2 * s * outer_sum_F1;
    *F2 = 1e-4 * s * s * outer_sum_F2;
}


static double
Iqabc(double qa, double qb, double qc,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double t)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;

    //HERE: Octahedron formula
    const double qx = qa * length_a;
    const double qy = qb * length_b;
    const double qz = qc * length_c;
    const double pf = 6./(qx*qx-qy*qy)*(1./(1.-(1.-t)*(1.-t)*(1.-t)));
    const double A = -qx*sin(qx)+0.5*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t));
    const double B =  qy*sin(qy)-0.5*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));
    const double AA = A/((qx*qx)-(qz*qz));
    const double BB = B/((qy*qy)-(qz*qz));

       // normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
    const double AP = pf*(AA+BB);
    // Multiply by contrast and volume
    const double s = (sld-solvent_sld) *(4./3.)* (length_a * length_b * length_c)*(1.-(1.-t)*(1.-t)*(1.-t));

    // Convert from [1e-12 A-1] to [cm-1]
    return 1.0e-4 * square(s * AP);
}
